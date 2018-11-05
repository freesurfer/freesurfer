from scipy.signal import fftconvolve
from scipy.special import expit
import os
import numpy as np
import freesurfer.gems as gems
from .utilities import requireNumpyArray, ensureDims
eps = np.finfo(float).eps


# Sample function
def sample(probs):
    return np.floor(probs + np.random.uniform(0, 1, np.shape(probs)))

# Gibbs sample function, k is the number of samples. Right now I'm using scipy convolutions, but in future one might
# want to use Tensorflow and GPUs to make it faster.
# TODO make it parallel or find a better solution (now it's sampling using just one core)
def gibbs_sample(k, x, W_conv, W_flip, bh, n_filter, paddings, n_hidden_1, n_hidden_2, n_hidden_3, width, height, depth):
    # Pre-allocate some of the values needed
    hk = np.zeros([n_hidden_1, n_hidden_2, n_hidden_3, n_filter])
    xk = np.zeros([width, height, depth, n_filter])
    for i in range(k):
        for f in range(n_filter):
            hk[:, :, :, f] = fftconvolve(x, W_conv[:, :, :, f], 'valid')
        hk = sample(expit(hk + bh))
        hk = np.pad(hk, paddings, 'constant')
        for f in range(n_filter):
            xk[:, :, :, f] = fftconvolve(hk[:, :, :, f], W_flip[:, :, :, f], 'valid')
        x = expit(np.sum(xk, 3))
    return x


# Function called by samsegment to model lesion shapes return FreeSurferLabels and volumesInCubicMm
def sampleCRBM(
        sampling_step,
        imageSize,
        modelPath,
        numberOfConstrast,
        wm_means,
        maskIndices,
        atlas,
        lesion_idx,
        dataMask,
        lesionsInit,
        likelihoods,
        savePath,
        imageToWorldTransformMatrix,
        FreeSurferLabels,
        nonCroppedImageSize,
        croppingOffset,
        lesionMask,
        names
):

    # Load the model
    model = np.load(modelPath)

    # Define the parameters of the model
    n_filter = model['numberOfFilter']
    conv_size = model['filterSize']

    width = imageSize[0]
    height = imageSize[1]
    depth = imageSize[2]

    n_hidden_1 = width - conv_size + 1
    n_hidden_2 = height - conv_size + 1
    n_hidden_3 = depth - conv_size + 1

    W_flip = model['weights']
    W_conv = W_flip[::-1, ::-1, ::-1, :]
    bh = model['hidBias']

    a1 = width - n_hidden_1
    a2 = height - n_hidden_2
    a3 = depth - n_hidden_3

    paddings = np.array([[a1, a1], [a2, a2], [a3, a3], [0, 0]], np.int8)

    # Create lesion and outlier masks
    # For mask if we have -1 mask below WM mean, +1 above, 0 nothing
    lesionsMask = np.zeros(imageSize, np.uint8)
    lesionsMask[lesionsInit] = 1
    mask = np.ones(imageSize, dtype=bool)
    k = 0
    for i in lesionMask:
        data = np.zeros(imageSize[0] * imageSize[1] * imageSize[2])
        data[np.reshape(maskIndices, [-1]) == 1] = dataMask[:, k]
        data = np.reshape(data, imageSize)
        if i == '-1':
            tmp = data < wm_means[k]
        elif i == '1':
            tmp = data > wm_means[k]
        elif i == '0':
            tmp = np.ones(imageSize, dtype=bool)
        mask = np.logical_and(mask, tmp)
        k = k + 1

    lesions = lesionsMask * mask

    lesion = np.nan_to_num(lesions, copy=False)

    posteriors_collected = np.zeros(atlas.shape, np.float32)
    non_lesion = np.where(np.arange(0, atlas.shape[1]) != lesion_idx)[0]

    for i in range(sampling_step):
        # Run network and get prior
        priors = atlas / 65535
        bv = priors[:, lesion_idx]
        data = gibbs_sample(1, lesion, W_conv, W_flip, bh, n_filter, paddings, n_hidden_1, n_hidden_2, n_hidden_3,
                            width, height, depth)

        # Assign lesion prior
        priors[:, lesion_idx] = np.reshape((data*mask)[maskIndices == 1], -1) * bv

        for l in non_lesion:  # don't loop over lesion
            priors[:, l] = priors[:, l] * (1 - priors[:, lesion_idx])

        # Normalize priors
        priors = np.nan_to_num(priors, copy=False)
        priors /= ensureDims(np.sum(priors, 1) + eps, 2)

        # Compute posteriors
        posteriors = np.multiply(priors, likelihoods)

        # Normalize posteriors
        posteriors = np.nan_to_num(posteriors, copy=False)
        posteriors /= ensureDims(np.sum(posteriors, 1) + eps, 2)

        # Multinomial sampling for posteriors
        print('Multinomial sampling')
        diff = np.cumsum(posteriors, axis=1) - np.random.random([posteriors.shape[0], 1])
        diff[diff < 0] = 100
        lesion_sample = np.argmin(diff, axis=1) == lesion_idx
        lesion_s = np.zeros(imageSize[0] * imageSize[1] * imageSize[2])
        lesion_s[np.reshape(maskIndices == 1, -1)] = lesion_sample
        lesion_s = np.reshape(lesion_s, imageSize)

        # Collect data
        posteriors_collected += posteriors

        # Prepare data to feed to the network for next loop cycle
        lesion = np.nan_to_num(lesion_s, copy=False)

        print('sample ' + str(i + 1) + ' times')

    # Perform majority voting and save images
    print('Perform majority voting')
    posteriors = posteriors_collected / sampling_step

    # Compute volumes in mm^3
    volumeOfOneVoxel = np.abs(np.linalg.det(imageToWorldTransformMatrix[0:3, 0:3]))
    volumesInCubicMm = (np.sum(posteriors, axis=0)) * volumeOfOneVoxel

    # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
    structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)
    freeSurferSegmentation = np.zeros(imageSize, dtype=np.uint16)
    FreeSurferLabels = np.array(FreeSurferLabels, dtype=np.uint16)
    freeSurferSegmentation[maskIndices == 1] = FreeSurferLabels[structureNumbers]

    # Write to file, remembering to un-crop the segmentation to the original image size
    uncroppedFreeSurferSegmentation = np.zeros(nonCroppedImageSize, dtype=np.float32)
    uncroppedFreeSurferSegmentation[croppingOffset[0]: imageSize[0] + croppingOffset[0],
    croppingOffset[1]: imageSize[1] + croppingOffset[1],
    croppingOffset[2]: imageSize[2] + croppingOffset[2]] = freeSurferSegmentation
    print('Writing out freesurfer segmentation')
    gems.KvlImage(requireNumpyArray(uncroppedFreeSurferSegmentation)).write(
        os.path.join(savePath, 'FinalSegmentation.nii'),
        gems.KvlTransform(requireNumpyArray(imageToWorldTransformMatrix))
    )

    # Not used now, but it might be useful to have all the probability maps for each class
    # Also write the posteriors
    #for classNumber in range(atlas.shape[1]):
    #    tmp = np.zeros(nonCroppedImageSize)
    #    tmp_post = np.zeros(imageSize)
    #    tmp_post[maskIndices == 1] = posteriors[:, classNumber]
    #    tmp[croppingOffset[0]: imageSize[0] + croppingOffset[0],
    #    croppingOffset[1]: imageSize[1] + croppingOffset[1],
    #    croppingOffset[2]: imageSize[2] + croppingOffset[2]] = tmp_post
    #    gems.KvlImage(requireNumpyArray(tmp)).write(
    #        os.path.join(savePath, str(names[classNumber]) + '_posterior.nii'), gems.KvlTransform(
    #            requireNumpyArray(imageToWorldTransformMatrix)))

    return [FreeSurferLabels, volumesInCubicMm]
