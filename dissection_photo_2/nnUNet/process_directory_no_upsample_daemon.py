import os

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import glob
import multiprocessing
import sys
import threading
import traceback
from argparse import ArgumentParser
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from time import sleep

import cv2
import numpy as np
import torch
sys.path.insert(0, os.path.dirname(__file__))
import nnunetv2.inference.predict_from_raw_data as nnunet_predict
from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor

EXIT_CALL = "ALL_DONE.txt"
POLL_INTERVAL = 0.2
MAX_ATTEMPTS = 2
MAX_POOL_RESTARTS = 3
STALE = "stale"
predictor = None

def log_error(message):
    print(message, file=sys.stderr, flush=True)

def parent_alive(parent_pid):
    if os.getppid() != parent_pid:
        return False
    try:
        os.kill(parent_pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    return True

def start_parent_watchdog(parent_pid):
    def watch():
        while parent_alive(parent_pid):
            sleep(1)
        try:
            log_error("Parent process is gone, exiting")
        finally:
            os._exit(1)

    threading.Thread(target=watch, daemon=True).start()

def case_name(image_file):
    return os.path.splitext(os.path.basename(image_file))[0]

def remove_file(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass

def image_stamp(image_file):
    try:
        st = os.stat(image_file)
    except OSError:
        return None
    return st.st_mtime_ns, st.st_ino

def close_call(image_file, out_dir, stamp):
    remove_file(image_file + ".txt")
    if image_stamp(image_file) == stamp:
        return True
    # the GUI re-queued the photo meanwhile: restore its call and drop the old verdict
    remove_file(os.path.join(out_dir, case_name(image_file) + ".failed"))
    open(image_file + ".txt", "a").close()
    return False

def mark_failed(image_file, out_dir, reason, stamp):
    with open(os.path.join(out_dir, case_name(image_file) + ".failed"), "w") as f:
        f.write(reason)
    kept = close_call(image_file, out_dir, stamp)
    if kept:
        log_error(f"Segmentation failed for {image_file}:\n{reason}")
    return kept

def load_downsampled(input_file):
    I = cv2.imread(input_file)
    if I is None:
        raise ValueError(f"Failed to read image from {input_file}")

    Id = cv2.resize(I, [0,0], fx = 0.2, fy = 0.2, interpolation=cv2.INTER_AREA)
    return np.stack([Id[:, :, 2], Id[:, :, 1], Id[:, :, 0]])[:, None]

def parse_calls(in_dir, running):
    calls = [call[:-len(".txt")] for call in glob.glob(os.path.join(glob.escape(in_dir), '*.txt'))
             if os.path.basename(call) != EXIT_CALL]
    return sorted(image_file for image_file in calls if image_file not in running)

def pick_device(preference):
    """auto: an NVIDIA GPU, else the Apple GPU, else the CPU. gpu: the same without
    the CPU, which is still used (with a message) when there is no GPU. cpu: the CPU."""
    if preference == "cpu":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    if preference == "gpu":
        log_error("No graphics card is available (NVIDIA or Apple GPU), using the CPU")
    return "cpu"

def folds_for(folds, device):
    """auto: all five folds on a GPU (CUDA or Apple), fold 0 on the CPU."""
    if folds == "5" or (folds == "auto" and device in ("cuda", "mps")):
        return (0, 1, 2, 3, 4)
    return (0,)

def build_predictor(model_path, device, folds):
    predictor = nnUNetPredictor(device=torch.device(device), allow_tqdm=False)
    predictor.initialize_from_trained_model_folder(model_path, use_folds=folds_for(folds, device))
    return predictor

def predict(image):
    # Inside inference mode: nnUNet sums the folds in place, which fails on the CPU
    # outside it ("Inplace update to inference tensor outside InferenceMode").
    with torch.inference_mode():
        _, probabilities = predictor.predict_single_npy_array(
            image, {'spacing': (999, 1, 1)}, save_or_return_probabilities=True)
    if not np.isfinite(probabilities).all():
        raise RuntimeError("The network produced non-finite values")
    return probabilities

def load_model(args):
    global predictor
    torch.set_num_threads(args.threads)
    # nnUNet caps prediction at its own default_num_processes (8): use ours.
    nnunet_predict.default_num_processes = args.threads

    device = pick_device(args.device)
    predictor = build_predictor(args.model_path, device, args.folds)
    if device == "mps":
        try:
            predict(np.zeros((3, 1, *predictor.configuration_manager.patch_size), np.float32))
        except Exception:
            log_error(f"The Apple GPU failed at startup, using the CPU:\n{traceback.format_exc()}")
            device = "cpu"
            predictor = build_predictor(args.model_path, device, args.folds)
    print(f"Segmenting on {device} with {len(folds_for(args.folds, device))} fold(s) and {args.threads} thread(s)", flush=True)

def segment(image_file, out_dir, stamp):
    probabilities = predict(load_downsampled(image_file))

    npz = os.path.join(out_dir, case_name(image_file) + ".npz")
    part = npz + ".part"
    try:
        with open(part, "wb") as f:
            np.savez(f, probabilities=probabilities)
        os.replace(part, npz)
    finally:
        remove_file(part)
    # only one job per image runs at a time, so this never removes a newer result
    if image_stamp(image_file) != stamp:
        remove_file(npz)
        return STALE

def segment_job(image_file, out_dir, stamp):
    try:
        return segment(image_file, out_dir, stamp)
    except Exception:
        return traceback.format_exc()

def init_worker(args, daemon_pid):
    start_parent_watchdog(daemon_pid)
    cv2.setNumThreads(1)
    load_model(args)

def start_workers(args):
    if args.workers == 1:
        load_model(args)
        return ThreadPoolExecutor(max_workers=1)

    for var in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = str(args.threads)

    executor = ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=init_worker,
        initargs=(args, os.getpid())
    )
    executor.submit(int)
    return executor

def main(args):
    start_parent_watchdog(os.getppid())
    in_dir, out_dir = args.in_dir, args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    try:
        executor = start_workers(args)
    except Exception:
        log_error(f"Could not load the model:\n{traceback.format_exc()}")
        return 1

    in_flight = {}
    stamps = {}
    attempts = {}
    solo = set()
    failures = restarts = 0

    def handle_result(image_file, error):
        nonlocal failures, restarts
        stamp = stamps[image_file]
        # the GUI re-wrote the image while it was being segmented: forget this run, the call is re-dispatched
        if error == STALE or image_stamp(image_file) != stamp:
            print(f"{image_file} changed while it was segmented, segmenting it again", flush=True)
        elif error is None:
            close_call(image_file, out_dir, stamp)
            restarts = 0
        elif attempts.get(image_file, 0) + 1 < MAX_ATTEMPTS:
            attempts[image_file] = attempts.get(image_file, 0) + 1
            log_error(f"Attempt {attempts[image_file]} failed for {image_file}, retrying:\n{error}")
            return
        elif mark_failed(image_file, out_dir, error, stamp):
            failures += 1
        attempts.pop(image_file, None)
        solo.discard(image_file)

    while True:
        finish_dir = os.path.isfile(os.path.join(in_dir, EXIT_CALL))
        pending = parse_calls(in_dir, set(in_flight.values()))
        if finish_dir and not pending and not in_flight:
            print("Finished processing images...")
            remove_file(os.path.join(in_dir, EXIT_CALL))
            return 1 if failures else 0

        crashed = False
        for image_file in pending:
            limit = 1 if image_file in solo or not solo.isdisjoint(in_flight.values()) else args.workers
            if len(in_flight) >= limit:
                break
            stamps[image_file] = image_stamp(image_file)
            try:
                in_flight[executor.submit(segment_job, image_file, out_dir, stamps[image_file])] = image_file
            except BrokenProcessPool:
                crashed = True

        if in_flight:
            done, _ = wait(in_flight, timeout=POLL_INTERVAL, return_when=FIRST_COMPLETED)
            for future in done:
                try:
                    error = future.result()
                except BrokenProcessPool:
                    crashed = True
                    continue
                handle_result(in_flight.pop(future), error)
        elif not crashed:
            sleep(POLL_INTERVAL)

        if crashed:
            solo.update(in_flight.values())
            for image_file in in_flight.values():
                handle_result(image_file, "A worker process stopped unexpectedly")
            in_flight.clear()
            restarts += 1
            if restarts > MAX_POOL_RESTARTS:
                reason = f"The worker processes stopped unexpectedly {restarts} times"
                for image_file in parse_calls(in_dir, set()):
                    mark_failed(image_file, out_dir, reason, image_stamp(image_file))
                return 1
            executor = start_workers(args)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--in_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--workers", type=int, default=1,
                        help="photos segmented at once, each by a process with its own copy of the model")
    parser.add_argument("--device", choices=("auto", "gpu", "cpu"), default="auto",
                        help="auto: a GPU if there is one, else the CPU; gpu: prefer the GPU; cpu: the CPU only")
    parser.add_argument("--folds", choices=("auto", "1", "5"), default="auto",
                        help="model folds to average: 1 (fold 0) or 5; auto: 5 on a GPU, 1 on the CPU")
    parser.add_argument("--max_threads", type=int, default=8,
                        help="CPU threads to use at most, shared by the workers; capped at what is available")

    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.max_threads < 1:
        parser.error("--max_threads must be at least 1")
    available = torch.get_num_threads()
    if args.max_threads > available:
        log_error(f"Asked for {args.max_threads} threads, but only {available} are available: using {available}")
    args.threads = max(1, min(args.max_threads, available) // args.workers)

    sys.exit(main(args))
