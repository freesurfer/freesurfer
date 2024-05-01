import os

import cv2
import numpy as np


def registration(
    true_width,
    true_height,
    template,
    des_template,
    centers,
    kp_template,
    input_image,
    output_dir=None,
    horizontal_ruler=None,
    vertical_ruler=None,
):
    # Constants
    DEBUG = False
    sift_res = 1024  # resolution at which SIFT operates;  TODO: make consistent with that of main.py
    reference_pixel_size = (
        0.1  # resolution of the final, perspective corrected image (in mm)
    )

    input_path, input_ext = os.path.splitext(input_image)
    _, input_name = os.path.split(input_path)
    output_image = input_name + "_corrected" + input_ext

    if output_dir is None:
        output_dir = os.getcwd()

    output_image = os.path.join(output_dir, output_image)

    # Read in new image to process
    target = cv2.imread(input_image, cv2.IMREAD_GRAYSCALE)
    target_fullsize_rgb = cv2.imread(input_image)

    # Resize image so smaller dimension is sift_res
    factor = sift_res / np.min(target.shape)
    new_target_size = np.round(np.flip(target.shape) * factor).astype(int)
    target = cv2.resize(target, dsize=new_target_size, interpolation=cv2.INTER_AREA)

    # Detect keypoints with SIFT
    sift = cv2.SIFT_create()
    kp_target, des_target = sift.detectAndCompute(target, None)

    if DEBUG:

        import matplotlib.pyplot as plt

        kp_im_template = template.copy()
        kp_im_template = cv2.drawKeypoints(
            template,
            kp_template,
            kp_im_template,
            flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS,
        )
        kp_im_target = target.copy()
        kp_im_target = cv2.drawKeypoints(
            target,
            kp_target,
            kp_im_target,
            flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS,
        )

        plt.figure(), plt.imshow(kp_im_template, aspect="equal"), plt.title(
            "Key points in template image"
        ), plt.show(block=False)
        plt.figure(), plt.imshow(kp_im_target, aspect="equal"), plt.title(
            "Key points in target image"
        ), plt.show(block=False)

    # Keypoint Matching
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)  # Brute force is fine

    # Match and extract points
    matches = bf.match(des_template, des_target)
    matches = sorted(matches, key=lambda x: x.distance)
    template_pts = np.float32([kp_template[m.queryIdx].pt for m in matches]).reshape(
        -1, 1, 2
    )
    target_pts = np.float32([kp_target[m.trainIdx].pt for m in matches]).reshape(
        -1, 1, 2
    )

    # Fit transform and apply to corner
    M, _ = cv2.findHomography(template_pts, target_pts, cv2.RANSAC, 2.0)
    centers_target = cv2.perspectiveTransform(centers.reshape(-1, 1, 2), M)

    if DEBUG:
        img = cv2.drawMatches(
            template,
            kp_template,
            target,
            kp_target,
            matches,
            None,
            flags=cv2.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS,
        )
        plt.figure(), plt.imshow(img, aspect="equal"), plt.title(
            "Matching key points"
        ), plt.show(block=False)
        img = cv2.polylines(
            target, [np.int32(centers_target)], True, 55, 3, cv2.LINE_AA
        )
        plt.figure(), plt.imshow(img, aspect="equal"), plt.title(
            "Detected corners"
        ), plt.show(block=False)

    # Now that we have detected the centers of the corners, we go back to the original coordinates
    centers_target = centers_target / factor
    target = target_fullsize_rgb

    # Now we only have to compute the final transform. The only caveat is the ordering of the corners...
    # We reorder then to NW, NE, SW, SE
    centers_target_reordered = np.zeros_like(centers_target)

    cost = centers_target[:, 0, 0] + centers_target[:, 0, 1]
    idx = np.argmin(cost)
    centers_target_reordered[0, 0, :] = centers_target[idx, 0, :]
    centers_target[idx, 0, :] = 0

    cost = -centers_target[:, 0, 0] + centers_target[:, 0, 1]
    cost[cost == 0] = 1e10
    idx = np.argmin(cost)
    centers_target_reordered[1, 0, :] = centers_target[idx, 0, :]
    centers_target[idx, 0, :] = 0

    cost = centers_target[:, 0, 0] - centers_target[:, 0, 1]
    cost[cost == 0] = 1e10
    idx = np.argmin(cost)
    centers_target_reordered[2, 0, :] = centers_target[idx, 0, :]
    centers_target[idx, 0, :] = 0

    cost = -centers_target[:, 0, 0] - centers_target[:, 0, 1]
    cost[cost == 0] = 1e10
    idx = np.argmin(cost)
    centers_target_reordered[3, 0, :] = centers_target[idx, 0, :]
    centers_target[idx, 0, :] = 0

    # We now define the target coordinates using the reerence resolution
    ref_coords = np.zeros_like(centers_target)

    ref_coords[0, 0, 0] = 0
    ref_coords[0, 0, 1] = 0

    ref_coords[1, 0, 0] = np.round(true_width / reference_pixel_size) - 1
    ref_coords[1, 0, 1] = 0

    ref_coords[2, 0, 0] = 0
    ref_coords[2, 0, 1] = np.round(true_height / reference_pixel_size) - 1

    ref_coords[3, 0, 0] = np.round(true_width / reference_pixel_size) - 1
    ref_coords[3, 0, 1] = np.round(true_height / reference_pixel_size) - 1

    PAD = 10.0 / reference_pixel_size  # pad 10mm (in pixels)
    ref_coords = ref_coords + PAD

    # We compute the final perspective transform
    M2, _ = cv2.findHomography(centers_target_reordered, ref_coords)
    deformed_image = cv2.warpPerspective(
        target,
        M2,
        (
            (PAD + ref_coords[1, 0, 0]).astype(int) + 1,
            (PAD + ref_coords[2, 0, 1]).astype(int) + 1,
        ),
    )

    # Add rulers if needed
    if horizontal_ruler is not None:
        image_with_ruler = np.zeros(
            (
                deformed_image.shape[0] + horizontal_ruler.shape[0],
                deformed_image.shape[1] + vertical_ruler.shape[1],
                3,
            ),
            dtype="uint8",
        )

        image_with_ruler[
            0 : deformed_image.shape[0], 0 : deformed_image.shape[1], :
        ] = deformed_image
        image_with_ruler[
            deformed_image.shape[0] :, 0 : -vertical_ruler.shape[1], :
        ] = horizontal_ruler[:, 0 : deformed_image.shape[1], :]
        image_with_ruler[
            0 : deformed_image.shape[0], -vertical_ruler.shape[1] :, :
        ] = vertical_ruler[0 : deformed_image.shape[0], :, :]

        deformed_image = image_with_ruler

    cv2.imwrite(output_image, deformed_image)

    if DEBUG:
        plt.figure(), plt.imshow(deformed_image, aspect="equal"), plt.title(
            "Perspective / pixel size corrected image"
        ), plt.show(block=False)



#### Function added by Eugenio
def fiducial_detection(
    true_width,
    true_height,
    template,
    des_template,
    centers,
    kp_template,
    input_image,
    output_file
):
    # Constants
    sift_res = 1024  # resolution at which SIFT operates;  TODO: make consistent with that of main.py
    reference_pixel_size = (
        0.1  # resolution of the final, perspective corrected image (in mm)
    )

    # Read in new image to process
    target = cv2.imread(input_image, cv2.IMREAD_GRAYSCALE)
    target_fullsize_rgb = cv2.imread(input_image)

    # Resize image so smaller dimension is sift_res
    factor = sift_res / np.min(target.shape)
    new_target_size = np.round(np.flip(target.shape) * factor).astype(int)
    target = cv2.resize(target, dsize=new_target_size, interpolation=cv2.INTER_AREA)

    # Detect keypoints with SIFT
    sift = cv2.SIFT_create()
    kp_target, des_target = sift.detectAndCompute(target, None)

    # Keypoint Matching
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)  # Brute force is fine

    # Match and extract points
    matches = bf.match(des_template, des_target)
    matches = sorted(matches, key=lambda x: x.distance)
    template_pts = np.float32([kp_template[m.queryIdx].pt for m in matches]).reshape(
        -1, 1, 2
    )
    target_pts = np.float32([kp_target[m.trainIdx].pt for m in matches]).reshape(
        -1, 1, 2
    )

    # Fit transform and apply to corner
    M, _ = cv2.findHomography(template_pts, target_pts, cv2.RANSAC, 2.0)
    centers_target = cv2.perspectiveTransform(centers.reshape(-1, 1, 2), M)

    # Now that we have detected the centers of the corners, we go back to the original coordinates
    centers_target = centers_target / factor

    # Write output to disk
    f = open(output_file, "w")
    f.write("TrueWidth: " + str(true_width) + '\n')
    f.write("TrueHeight: " + str(true_height) + '\n')
    for c in range(4):
        f.write('Corner_' + str(c+1) + ': ' + str(centers_target[c,0,0]) + ',' + str(centers_target[c,0,1]) + '\n')
    f.close()
