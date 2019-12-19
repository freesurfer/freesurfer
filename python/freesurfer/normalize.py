'''
These functions are rather specific to Amod's deeplearning code
and should probably be either moved or completely deleted. It might
make sense to develop some more universal normalization utility functions.
'''

import numpy as np


def piecewise_linear_normalize(in_img_data, ref_img_data):
    '''Function to piecewise linearly scale image intensities to training data landmarks.'''
    import sklearn.mixture
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
    print(ref_landmarks)
    print(in_landmarks)
    out_img_data = np.zeros(in_img_data.shape)

    # map intensities using these landmarks
    for i in range(len(in_landmarks)-1):
        m = (ref_landmarks[i+1] - ref_landmarks[i])/(in_landmarks[i+1] - in_landmarks[i])
        c = (in_landmarks[i+1]*ref_landmarks[i] - in_landmarks[i]*ref_landmarks[i+1])/(in_landmarks[i+1] - in_landmarks[i])
        out_img_data[(in_img_data > in_landmarks[i]) & (in_img_data <= in_landmarks[i+1])] = \
            m * in_img_data[(in_img_data > in_landmarks[i]) & (in_img_data <= in_landmarks[i+1])] + c

    out_img_data[(in_img_data > in_landmarks[-1])] = 255
    return out_img_data


def wm_peak_normalize(in_img_data):
    '''Function to scale image intensities by setting wm peak to 200.'''
    import sklearn.mixture
    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    # clf = mixture.GMM(n_components=3, covariance_type='full')
    clf = mixture.GaussianMixture(n_components=3, covariance_type='full', n_init=5)
    clf.fit(in_img_fg)
    # max of means is the wm centroid for t1w images
    wm_peak_intensity  = clf.means_.max()
    wm_scaling = 200.0 / wm_peak_intensity
    print(wm_peak_intensity)
    out_img_data = in_img_data * wm_scaling
    return out_img_data


def robust_normalize(in_img_data):
    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p01 = np.percentile(in_img_fg, q=1)
    p999 = np.percentile(in_img_fg, q=99)
    scaling = 255.0 / (p999 - p01)  # set p99 to 255
    in_img_data[(in_img_data < p01) & (in_img_data > 0 )] = p01
    print(scaling)
    out_img_data = (in_img_data) * scaling
    return out_img_data

def wm_peak_normalize_t2w(in_img_data):
    '''Function to scale image intensities by setting wm peak to 200.'''
    import sklearn.mixture
    in_img_flat = np.ravel(in_img_data, 'C')

    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p95 = np.percentile(in_img_fg, q=90)
    p05 = np.percentile(in_img_fg, q=10)

    in_img_fg = in_img_fg[(in_img_fg < p95) & (in_img_fg > p05)]
    in_img_fg = in_img_fg.reshape(-1,1)

    # clf = mixture.GMM(n_components=3, covariance_type='full')
    clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(in_img_fg)
    print('GMM centroids are ')
    print(sorted(clf.means_))
    wm_peak_intensity = sorted(clf.means_)[0]

    # h, bin_edges = np.histogram(in_img_fg, 500)
    # max_bin = np.argmax(h)
    # mode_h = max_bin * (bin_edges[1] - bin_edges[0])

    # max of means is the wm centroid for t1w images
    # wm_peak_intensity = mode_h
    wm_scaling = 0.3 / wm_peak_intensity
    print(wm_peak_intensity)

    out_img_data = in_img_data * wm_scaling
    return out_img_data

def wm_peak_normalize_t2w(in_img_data):
    """Function to scale image intensities by setting wm peak to 200"""
    import sklearn.mixture

    in_img_flat = np.ravel(in_img_data, 'C')

    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p95 = np.percentile(in_img_fg, q=90)
    p05 = np.percentile(in_img_fg, q=10)

    in_img_fg = in_img_fg[(in_img_fg < p95) & (in_img_fg > p05)]
    in_img_fg = in_img_fg.reshape(-1,1)

    # clf = mixture.GMM(n_components=3, covariance_type='full')
    clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(in_img_fg)
    print('GMM centroids are:')
    print(sorted(clf.means_))
    wm_peak_intensity = sorted(clf.means_)[0]
    #
    #
    # h, bin_edges = np.histogram(in_img_fg, 500)
    # max_bin = np.argmax(h)
    # mode_h = max_bin * (bin_edges[1] - bin_edges[0])

    # max of means is the wm centroid for t1w images
    #     wm_peak_intensity  = mode_h
    wm_scaling = 0.3 / wm_peak_intensity
    print(wm_peak_intensity)

    out_img_data = in_img_data * wm_scaling
    return out_img_data


def max_normalize(in_img_data):
    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p01 = np.percentile(in_img_fg, q=1)
    p999 = np.percentile(in_img_fg, q=99)
    scaling = 255.0 / (p999 - p01)  # set p99 to 255
    in_img_data[in_img_data < p01] = p01
    print(scaling)
    out_img_data = (in_img_data) * scaling
    return out_img_data


def histmatch(in_img_data, ref_img_data):
    # in_img_data = wm_peak_normalize(in_img_data)
    # ref_img_data = wm_peak_normalize(ref_img_data)
    in_img_data_flat = in_img_data.flatten()
    in_img_fg = in_img_data_flat[in_img_data_flat > 0]  # foreground is > 0

    ref_img_data_flat = ref_img_data.flatten()
    ref_img_fg = ref_img_data_flat[ref_img_data_flat > 0]  # foreground is > 0

    bins_in = np.linspace(0, 1, 255 / 1)
    bins_ref = np.linspace(0, 1, 255 / 1)

    hist_in = np.histogram(in_img_fg, bins=bins_in, range=(bins_in.min(), bins_in.max()))
    n_in = hist_in[0]
    bins_in = hist_in[1]

    hist_ref = np.histogram(ref_img_fg, bins=bins_ref, range=(bins_ref.min(), bins_ref.max()))
    n_ref = hist_ref[0]
    bins_ref = hist_ref[1]

    cdf_in_img = np.float64(np.cumsum(n_in))
    cdf_in_img = np.divide(cdf_in_img, cdf_in_img[-1])

    cdf_ref_img = np.float64(np.cumsum(n_ref))
    cdf_ref_img = np.divide(cdf_ref_img, cdf_ref_img[-1])

    interp_ref_values = np.interp(cdf_in_img, cdf_ref_img, bins_ref[1:])
    bins_in_z = np.append(0, bins_in)

    out_img_data = np.copy(in_img_data)
    for i in range(1, len(bins_in)):
        out_img_data[(in_img_data > bins_in_z[i - 1]) & (in_img_data <= bins_in_z[i])] = interp_ref_values[i - 1]

    return out_img_data, bins_in_z, interp_ref_values
