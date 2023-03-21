import argparse
import glob
import os
import sys

import cv2
import numpy as np
from registration import registration


def fiducials_correction(args):
    """Perform correction of photos with fiducials"""

    npz_file_path = args.npz_file
    input_folder_path = args.in_dir
    output_folder_path = args.out_dir

    if not os.path.isfile(args.npz_file):
        raise Exception("Warning", "Calibration File not Found")

    if not os.path.isdir(input_folder_path):
        raise Exception("Warning", "Input Directory Cannot be Empty")

    if not os.path.isdir(output_folder_path):
        os.makedirs(output_folder_path, exist_ok=True)

    # Read data from model file
    variables = np.load(npz_file_path, allow_pickle=True)
    true_width = variables["true_w"].astype("float")
    true_height = variables["true_h"].astype("float")
    template = variables["img_template"]
    kp_template_tmp = variables["kp_template"]
    des_template = variables["des_template"]
    centers = variables["centers"]

    # reassemble the key points (which we split for saving to disk)
    kp_template = []
    for point in kp_template_tmp:
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

    input_images = sorted(glob.glob(os.path.join(input_folder_path, "*.*")))

    horizontal_ruler = cv2.imread("./resources/horizontal.png")
    vertical_ruler = cv2.imread("./resources/vertical.png")

    for input_image in input_images:
        try:
            registration(
                true_width,
                true_height,
                template,
                des_template,
                centers,
                kp_template,
                input_image,
                output_folder_path,
                horizontal_ruler,
                vertical_ruler,
            )
        except:
            print(f"failed on {input_image}")

    print("Performed Registration Successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_dir", type=str, dest="in_dir", default=None)
    parser.add_argument("--calibration_file", type=str, dest="npz_file", default=None)
    parser.add_argument("--out_dir", type=str, dest="out_dir", default=None)

    # If running the code in debug mode
    gettrace = getattr(sys, "gettrace", None)

    if gettrace():
        sys.argv = [
            "func_fiducials_correction.py",
            "--in_dir",
            "/space/calico/1/users/Harsha/photo-calibration-gui/misc/fiducials_correction_input/",
            "--calibration_file",
            "/tmp/cal.npz",
            "--out_dir",
            "/tmp",
        ]

    args = parser.parse_args()

    fiducials_correction(args)
