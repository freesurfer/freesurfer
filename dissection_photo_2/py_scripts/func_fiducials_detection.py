import argparse
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
import cv2
import numpy as np
from registration import fiducial_detection
from func_retrospective_correction import retrospective_correction


def calibrated_correction(args):
    """Detect the fiducials in each photo and correct it with them.

    One run for all the photos: the camera may move between photos, so every
    photo gets its own detection. A few photos are worked on at a time (OpenCV
    and the image codecs release the GIL). Progress goes to stdout, one
    tab-separated line per event, for MainWindow.cpp (SIZE, START, DONE, FAILED).
    """

    if not os.path.isfile(args.npz_file):
        raise Exception("Warning", "Calibration File not Found")

    variables = np.load(args.npz_file, allow_pickle=True)
    true_width = float(variables["true_w"])
    true_height = float(variables["true_h"])
    template = variables["img_template"]
    des_template = variables["des_template"]
    centers = variables["centers"]

    # reassemble the key points (which we split for saving to disk)
    kp_template = []
    for point in variables["kp_template"]:
        temp = cv2.KeyPoint(
            x=point[0],
            y=point[1],
            size=point[2],
            angle=point[3],
            response=point[4],
            octave=int(point[5]),
            class_id=int(point[6]),
        )
        kp_template.append(temp)

    lock = threading.Lock()

    def report(*fields):
        with lock:  # whole lines, however the threads interleave
            sys.stdout.write("\t".join(fields) + "\n")
            sys.stdout.flush()

    def correct(image):
        report("START", image)
        try:
            corners = fiducial_detection(
                true_width, true_height, template, des_template, centers, kp_template, image
            )
            retrospective_correction(
                argparse.Namespace(
                    in_img=image,
                    pos_tuple=corners.tolist(),
                    e1=true_width,
                    e2=true_height,
                    out_dir=args.out_dir,
                )
            )
            report("DONE", image, " ".join(f"{v:.1f}" for v in corners.ravel()))
        except Exception as e:  # this photo is left for the user to correct by hand
            report("FAILED", image, str(e))

    report("SIZE", str(true_width), str(true_height))
    with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
        list(pool.map(correct, args.in_images))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_images", type=str, nargs="+", dest="in_images", default=None)
    parser.add_argument("--calibration_file", type=str, dest="npz_file", default=None)
    parser.add_argument("--out_dir", type=str, dest="out_dir", default=None)

    args = parser.parse_args()

    calibrated_correction(args)
