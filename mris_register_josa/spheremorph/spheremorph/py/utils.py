import numpy as np


# adapted from freesurfer deeplearn
def norm_curvature(curv, which_norm='Median', norm_percentile=97, std_thresh=3):
    if which_norm == 'Percentile':
        normed = np.clip(curv / np.percentile(curv, norm_percentile), 0, 1)
    elif which_norm == 'Median':
        min_clip = np.percentile(curv, 100 - norm_percentile)
        max_clip = np.percentile(curv, norm_percentile)
        st = np.std(np.clip(curv, min_clip, max_clip))
        normed = np.clip(((curv - np.median(curv)) / st), -std_thresh, std_thresh)
    else:
        normed = (curv - curv.mean()) / np.std(curv)

    return normed
