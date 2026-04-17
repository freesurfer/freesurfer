import glob
import os
import sys
from argparse import ArgumentParser

import numpy as np
from skimage.io import imread
from skimage.measure import label as bwlabel


def mask_to_cc(args):
    """Label connected components in 2-D binary image"""

    if not os.path.isdir(args.in_dir):
        print("Input directory does not exist")
        sys.exit()

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    mask_files = sorted(glob.glob(os.path.join(args.in_dir, "*")))
    if not mask_files:
        print("No masks found in the input directory")
        sys.exit()

    for current_mask in mask_files:
        # Open mask
        try:
            mask = imread(current_mask, as_gray=True)
        except:
            continue

        mask = mask > np.max(mask) / 2

        # Create connected components
        connected_components = bwlabel(mask)

        # Save connected components
        file_name = os.path.split(current_mask)[-1]
        file_name, _ = os.path.splitext(file_name)
        np.savez_compressed(os.path.join(args.out_dir, file_name), cc=connected_components.astype('uint16'))


if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--in_dir", type=str, dest="in_dir", default=None)
    parser.add_argument("--out_dir", type=str, dest="out_dir", default=None)

    # If running the code in debug mode
    gettrace = getattr(sys, "gettrace", None)

    if gettrace():
        sys.argv = [
            "func_mask_to_cc.py",
            "--in_dir",
            "/space/calico/1/users/Harsha/photo-calibration-gui/misc/masked",
            "--out_dir",
            "/space/calico/1/users/Harsha/photo-calibration-gui/misc/cc_temp/",
        ]

    args = parser.parse_args()

    mask_to_cc(args)

    # example call:
    # fspython func_mask_to_cc.py \
    #   --in_dir  "/space/calico/1/users/Harsha/photo-calibration-gui/misc/masked" \
    #   --out_dir  "/space/calico/1/users/Harsha/photo-calibration-gui/misc/cc_temp"
