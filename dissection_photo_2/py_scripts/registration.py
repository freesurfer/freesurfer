import cv2
import numpy as np


def fiducial_detection(
    true_width,
    true_height,
    template,
    des_template,
    centers,
    kp_template,
    input_image,
):
    """Centers of the four fiducials in input_image, in its own pixel coordinates (4 x 2)."""
    sift_res = 1024  # resolution at which SIFT operates; must match func_fiducials_calibration.py

    target = cv2.imread(input_image, cv2.IMREAD_GRAYSCALE)

    # Resize image so smaller dimension is sift_res
    factor = sift_res / np.min(target.shape)
    new_target_size = np.round(np.flip(target.shape) * factor).astype(int)
    target = cv2.resize(target, dsize=new_target_size, interpolation=cv2.INTER_AREA)

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
    # Any 4 matches fit a transform exactly, so a photo without the markers still
    # "fits": ask for clearly more matches that agree with it.
    min_inliers = 10
    M, inliers = (None, None) if len(matches) < 4 else cv2.findHomography(template_pts, target_pts, cv2.RANSAC, 2.0)
    if M is None or int(inliers.sum()) < min_inliers:
        raise ValueError("the calibration markers were not found")
    centers_target = cv2.perspectiveTransform(centers.reshape(-1, 1, 2), M)

    # Now that we have detected the centers of the corners, we go back to the original coordinates
    centers_target = centers_target / factor

    return centers_target[:, 0, :]
