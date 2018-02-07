import numpy as np
import nibabel as nib
import skimage
from sklearn.neighbors import KernelDensity
from sklearn import mixture


def piecewise_linear_normalize(in_img_data, ref_img_data):
    """Function to piecewise linearly scale image intensities to training data landmarks"""

    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    clf_in = mixture.GaussianMixture(n_components=3, covariance_type='full')
    clf_in.fit(in_img_fg)

    ref_img_flat = np.ravel(ref_img_data, 'C')
    ref_img_fg = ref_img_flat[ref_img_flat > 0].reshape(-1, 1)
    clf_ref = mixture.GaussianMixture(n_components=3, covariance_type='full')
    clf_ref.fit(ref_img_fg)

    in_landmarks = np.asarray(sorted(clf_in.means_.squeeze()))
    in_wm_std = np.sqrt(clf_in.covariances_[np.argmax(clf_in.means_)])
    in_wm_threshold = in_landmarks[2] + 2*in_wm_std[0]
    in_landmarks = np.append(np.asarray([0]), np.append(in_landmarks, in_wm_threshold))

    ref_landmarks = np.asanyarray(sorted(clf_ref.means_.squeeze()))
    ref_wm_std = np.sqrt(clf_in.covariances_[np.argmax(clf_in.means_)])
    ref_wm_threshold = 255
    ref_landmarks = np.append(np.asarray([0]), np.append(ref_landmarks, ref_wm_threshold))
    print ref_landmarks
    print in_landmarks
    out_img_data = np.zeros(in_img_data.shape)
    # map intensities using these landmarks
    for i in range(len(in_landmarks)-1):
        m = (ref_landmarks[i+1] - ref_landmarks[i])/(in_landmarks[i+1] - in_landmarks[i])
        c = (in_landmarks[i+1]*ref_landmarks[i] - in_landmarks[i]*ref_landmarks[i+1])/(in_landmarks[i+1] - in_landmarks[i])

        out_img_data[(in_img_data > in_landmarks[i]) & (in_img_data <= in_landmarks[i+1])] = \
            m*in_img_data[(in_img_data > in_landmarks[i]) & (in_img_data <= in_landmarks[i+1])] + c

    out_img_data[(in_img_data > in_landmarks[-1])] = 255
    return out_img_data



def wm_peak_normalize(in_img_data):
    """Function to scale image intensities by setting wm peak to 200"""


    in_img_flat = np.ravel(in_img_data, 'C')

    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    # clf = mixture.GMM(n_components=3, covariance_type='full')

    clf = mixture.GaussianMixture(n_components=3, covariance_type='full')
    clf.fit(in_img_fg)
    # max of means is the wm centroid for t1w images
    wm_peak_intensity  = clf.means_.max()
    wm_scaling = 200 / wm_peak_intensity
    print(wm_peak_intensity)

    out_img_data = in_img_data * wm_scaling
    return out_img_data


def robust_normalize(in_img_data):
    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p99 = np.percentile(in_img_fg, q=99)

    # set p99 to 255
    scaling = 255.0 / p99
    print(scaling)
    out_img_data = in_img_data * scaling
    return out_img_data



