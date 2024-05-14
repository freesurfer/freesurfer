import argparse
import os
import cv2
import numpy as np
from registration import fiducial_detection


def fiducials_detection(args):
    """Perform automated detection of fiducials in a photo"""

    npz_file_path = args.npz_file
    input_image = args.in_image
    output_file = args.out_file

    if not os.path.isfile(npz_file_path):
        raise Exception("Warning", "Calibration File not Found")

    if not os.path.isfile(input_image):
        raise Exception("Warning", "Input Directory Cannot be Empty")

    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

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

    fiducial_detection(
        true_width,
        true_height,
        template,
        des_template,
        centers,
        kp_template,
        input_image,
        output_file
    )

    print("Performed fiducial detection")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_image", type=str, dest="in_image", default=None)
    parser.add_argument("--calibration_file", type=str, dest="npz_file", default=None)
    parser.add_argument("--out_file", type=str, dest="out_file", default=None)

    args = parser.parse_args()

    fiducials_detection(args)