import argparse
import os
import sys

import cv2
import numpy as np
from PIL import Image, ImageOps
from scipy.optimize import linear_sum_assignment


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.split(values))

    def chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    def split(self, s):
        try:
            s = list(map(float, s))
        except:
            s = list(map(float, s[0].split()))

        coords = list(self.chunks(s, 2))
        if len(coords[-1]) != 2:
            print("Invalid coordinates")
            sys.exit()
        return coords


def retrospective_correction(args):
    """Warp in_img to 0.1 mm/pixel from 2 (scale bar), 3 or 4 (corner) clicks, with rulers."""

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)

    # MainWindow.cpp rewrites these two paths in the script text: keep them byte-identical
    args.horizontal_ruler = cv2.imread("./resources/horizontal.png")
    args.vertical_ruler = cv2.imread("./resources/vertical.png")

    input_path, input_ext = os.path.splitext(args.in_img)
    _, input_name = os.path.split(input_path)
    args.out_img = os.path.join(args.out_dir, input_name + "_deformed" + input_ext)

    args.img_fullres = Image.open(args.in_img)
    args.img_fullres = ImageOps.exif_transpose(args.img_fullres)

    n_points = len(args.pos_tuple)

    if n_points in [3, 4]:
        true_width, true_height = args.e1, args.e2
    elif n_points == 2:
        # We pretend the user clicked on the 4 corners of the image
        # and make the true_width and true_height proportional to the provided length
        pix_dist = np.sqrt(
            (args.pos_tuple[0][0] - args.pos_tuple[1][0]) ** 2
            + (args.pos_tuple[0][1] - args.pos_tuple[1][1]) ** 2
        )
        pix_siz = float(args.e1) / pix_dist

        true_width = pix_siz * args.img_fullres.size[0]
        true_height = pix_siz * args.img_fullres.size[1]

        args.pos_tuple = [
            [0, 0],
            [0, args.img_fullres.size[1] - 1],
            [args.img_fullres.size[0] - 1, 0],
            [args.img_fullres.size[0] - 1, args.img_fullres.size[1] - 1],
        ]
    else:
        raise ValueError("Wrong number of clicks")

    reference_pixel_size = 0.1

    centers_target = np.array(args.pos_tuple, dtype=float)[:, np.newaxis, :]

    # Match the clicks to the corners NW, NE, SW, SE (Hungarian algorithm)
    costNW = centers_target[:, 0, 0] + centers_target[:, 0, 1]
    costNW -= np.min(costNW)
    costNE = -centers_target[:, 0, 0] + centers_target[:, 0, 1]
    costNE -= np.min(costNE)
    costSW = centers_target[:, 0, 0] - centers_target[:, 0, 1]
    costSW -= np.min(costSW)
    costSE = -centers_target[:, 0, 0] - centers_target[:, 0, 1]
    costSE -= np.min(costSE)

    C = np.stack([costNW, costNE, costSW, costSE], axis=-1)
    if n_points == 3:
        C = np.concatenate([C, -100000 * np.ones([1, 4])], axis=0)

    _, idx = linear_sum_assignment(C)

    # Target corners NW, NE, SW, SE at the reference resolution, padded by 10 mm
    w = np.round(true_width / reference_pixel_size) - 1
    h = np.round(true_height / reference_pixel_size) - 1
    PAD = 10.0 / reference_pixel_size
    ref_coords = np.array([[[0, 0]], [[w, 0]], [[0, h]], [[w, h]]]) + PAD

    if n_points == 3:
        M2 = cv2.getAffineTransform(
            centers_target.astype(np.float32),
            ref_coords[idx[:n_points], :].astype(np.float32),
        )
        warp_function = cv2.warpAffine
    else:
        M2, _ = cv2.findHomography(centers_target, ref_coords[idx, :])
        warp_function = cv2.warpPerspective

    args.deformed_image = warp_function(
        np.asarray(args.img_fullres),
        M2,
        (
            (PAD + ref_coords[1, 0, 0]).astype(int) + 1,
            (PAD + ref_coords[2, 0, 1]).astype(int) + 1,
        ),
    )

    image_with_ruler = np.zeros(
        (
            args.deformed_image.shape[0] + args.horizontal_ruler.shape[0],
            args.deformed_image.shape[1] + args.vertical_ruler.shape[1],
            3,
        ),
        dtype="uint8",
    )

    image_with_ruler[
        0 : args.deformed_image.shape[0], 0 : args.deformed_image.shape[1], :
    ] = cv2.cvtColor(args.deformed_image, cv2.COLOR_RGB2BGR)
    image_with_ruler[
        args.deformed_image.shape[0] :, 0 : -args.vertical_ruler.shape[1], :
    ] = args.horizontal_ruler[:, 0 : args.deformed_image.shape[1], :]
    image_with_ruler[
        0 : args.deformed_image.shape[0], -args.vertical_ruler.shape[1] :, :
    ] = args.vertical_ruler[0 : args.deformed_image.shape[0], :, :]

    cv2.imwrite(args.out_img, image_with_ruler)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_img", type=str, dest="in_img", default=None)
    parser.add_argument("--points", nargs="+", dest="pos_tuple", action=SplitArgs)
    parser.add_argument("--width", nargs="?", type=float, dest="e1", default=None)
    parser.add_argument("--height", nargs="?", type=float, dest="e2", default=None)
    parser.add_argument("--out_dir", type=str, dest="out_dir", default=None)

    args = parser.parse_args()

    retrospective_correction(args)
