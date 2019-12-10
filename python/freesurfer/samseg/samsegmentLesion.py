from freesurfer.samseg.samsegment import samsegment, getModelSpecifications, getFullHyperparameters, \
    getLikelihoods, fitGMMParameters, getGaussianLikelihoods, evaluateMinLogPriorOfGMMParameters
from freesurfer.samseg.figures import initVisualizer
import freesurfer.gems as gems
from freesurfer.samseg.utilities import Specification
import numpy as np
import tensorflow as tf
from scipy.ndimage.interpolation import affine_transform
from scipy.stats import invwishart
from freesurfer.gems import kvlReadSharedGMMParameters
import os

# eps = np.spacing(1)
eps = np.finfo(float).eps


def getClassNumber(structureSearchString):
    #
    if structureSearchString is None:
        return None

    #
    global atlasDir_

    sharedGMMParameters = kvlReadSharedGMMParameters(os.path.join(atlasDir_, 'sharedGMMParameters.txt'))

    for classNumber, mergeOption in enumerate(sharedGMMParameters):
        for searchString in mergeOption.searchStrings:
            if structureSearchString in searchString:
                structureClassNumber = classNumber

    return structureClassNumber


# Encoder network
def get_encoder(data):
    epsilon = 1e-3

    graph = tf.get_default_graph()

    # First convolution
    conv_1_bias = graph.get_tensor_by_name('aconv_1/bias:0')
    conv_1_kernel = graph.get_tensor_by_name('aconv_1/kernel:0')
    conv_1 = tf.nn.conv3d(data, conv_1_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_1_bias
    conv_1 = tf.nn.relu(conv_1)
    mu_1 = graph.get_tensor_by_name('abatch_n_e_1/moving_mean:0')
    var_1 = graph.get_tensor_by_name('abatch_n_e_1/moving_variance:0')
    beta_1 = graph.get_tensor_by_name('abatch_n_e_1/beta:0')
    gamma_1 = graph.get_tensor_by_name('abatch_n_e_1/gamma:0')
    batch_n_e_1 = gamma_1 * ((conv_1 - mu_1) / tf.sqrt(var_1 + epsilon)) + beta_1

    # Second convolution
    conv_2_bias = graph.get_tensor_by_name('aconv_2/bias:0')
    conv_2_kernel = graph.get_tensor_by_name('aconv_2/kernel:0')
    conv_2 = tf.nn.conv3d(batch_n_e_1, conv_2_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_2_bias
    conv_2 = tf.nn.relu(conv_2)
    mu_2 = graph.get_tensor_by_name('abatch_n_e_2/moving_mean:0')
    var_2 = graph.get_tensor_by_name('abatch_n_e_2/moving_variance:0')
    beta_2 = graph.get_tensor_by_name('abatch_n_e_2/beta:0')
    gamma_2 = graph.get_tensor_by_name('abatch_n_e_2/gamma:0')
    batch_n_e_2 = gamma_2 * ((conv_2 - mu_2) / tf.sqrt(var_2 + epsilon)) + beta_2

    # Third convolution
    conv_3_bias = graph.get_tensor_by_name('aconv_3/bias:0')
    conv_3_kernel = graph.get_tensor_by_name('aconv_3/kernel:0')
    conv_3 = tf.nn.conv3d(batch_n_e_2, conv_3_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_3_bias
    conv_3 = tf.nn.relu(conv_3)
    mu_3 = graph.get_tensor_by_name('abatch_n_e_3/moving_mean:0')
    var_3 = graph.get_tensor_by_name('abatch_n_e_3/moving_variance:0')
    beta_3 = graph.get_tensor_by_name('abatch_n_e_3/beta:0')
    gamma_3 = graph.get_tensor_by_name('abatch_n_e_3/gamma:0')
    batch_n_e_3 = gamma_3 * ((conv_3 - mu_3) / tf.sqrt(var_3 + epsilon)) + beta_3

    # Fourth convolution
    conv_4_bias = graph.get_tensor_by_name('aconv_4/bias:0')
    conv_4_kernel = graph.get_tensor_by_name('aconv_4/kernel:0')
    conv_4 = tf.nn.conv3d(batch_n_e_3, conv_4_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_4_bias
    var_4 = graph.get_tensor_by_name('abatch_n_e_4/moving_variance:0')
    beta_4 = graph.get_tensor_by_name('abatch_n_e_4/beta:0')
    gamma_4 = graph.get_tensor_by_name('abatch_n_e_4/gamma:0')
    conv_4 = tf.nn.relu(conv_4)
    mu_4 = graph.get_tensor_by_name('abatch_n_e_4/moving_mean:0')
    batch_n_e_4 = gamma_4 * ((conv_4 - mu_4) / tf.sqrt(var_4 + epsilon)) + beta_4

    # Fifth convolution
    conv_5_bias = graph.get_tensor_by_name('aconv_5/bias:0')
    conv_5_kernel = graph.get_tensor_by_name('aconv_5/kernel:0')
    conv_5 = tf.nn.conv3d(batch_n_e_4, conv_5_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + conv_5_bias
    conv_5 = tf.nn.relu(conv_5)
    mu_5 = graph.get_tensor_by_name('abatch_n_e_5/moving_mean:0')
    var_5 = graph.get_tensor_by_name('abatch_n_e_5/moving_variance:0')
    beta_5 = graph.get_tensor_by_name('abatch_n_e_5/beta:0')
    gamma_5 = graph.get_tensor_by_name('abatch_n_e_5/gamma:0')
    batch_n_e_5 = gamma_5 * ((conv_5 - mu_5) / tf.sqrt(var_5 + epsilon)) + beta_5

    # Convolve to mu e sigma instead of using fully connected
    mu_bias = graph.get_tensor_by_name('amu/bias:0')
    mu_kernel = graph.get_tensor_by_name('amu/kernel:0')
    mu = tf.nn.conv3d(batch_n_e_5, mu_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + mu_bias

    sigma_bias = graph.get_tensor_by_name('asigma/bias:0')
    sigma_kernel = graph.get_tensor_by_name('asigma/kernel:0')
    sigma = tf.nn.conv3d(batch_n_e_5, sigma_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + sigma_bias
    sigma = tf.nn.softplus(sigma)

    return mu, sigma


# decoder network
def get_decoder(code, imageSize):
    epsilon = 1e-3

    graph = tf.get_default_graph()

    code_size = tf.shape(code)

    # First deconv layer
    deconv_0_bias = graph.get_tensor_by_name('adeconv_0/bias:0')
    deconv_0_kernel = graph.get_tensor_by_name('adeconv_0/kernel:0')
    # Deconvolution shape for VALID = stride * (input - 1) + kernel size
    deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3,
                             code_size[3] - 1 + 3, 16])

    hidden_0_dec = tf.nn.conv3d_transpose(code, deconv_0_kernel, output_shape=deconv_shape,
                                          strides=[1, 1, 1, 1, 1], padding='VALID') + deconv_0_bias

    mu_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_mean:0')
    var_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_variance:0')
    beta_1 = graph.get_tensor_by_name('abatch_n_d_0/beta:0')
    gamma_1 = graph.get_tensor_by_name('abatch_n_d_0/gamma:0')
    batch_n_d_1 = gamma_1 * ((hidden_0_dec - mu_1) / tf.sqrt(var_1 + epsilon)) + beta_1

    # Second deconv layer
    code_size = tf.shape(batch_n_d_1)
    deconv_1_bias = graph.get_tensor_by_name('adeconv_1/bias:0')
    deconv_1_kernel = graph.get_tensor_by_name('adeconv_1/kernel:0')
    deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3,
                             code_size[3] - 1 + 3, 16])

    hidden_2_dec = tf.nn.conv3d_transpose(batch_n_d_1, deconv_1_kernel, output_shape=deconv_shape,
                                          strides=[1, 1, 1, 1, 1], padding='VALID') + deconv_1_bias

    hidden_2_dec = tf.nn.relu(hidden_2_dec)
    mu_2 = graph.get_tensor_by_name('abatch_n_d_1/moving_mean:0')
    var_2 = graph.get_tensor_by_name('abatch_n_d_1/moving_variance:0')
    beta_2 = graph.get_tensor_by_name('abatch_n_d_1/beta:0')
    gamma_2 = graph.get_tensor_by_name('abatch_n_d_1/gamma:0')
    batch_n_d_2 = gamma_2 * ((hidden_2_dec - mu_2) / tf.sqrt(var_2 + epsilon)) + beta_2

    # Third deconv layer
    code_size = tf.shape(batch_n_d_2)
    deconv_2_bias = graph.get_tensor_by_name('adeconv_2/bias:0')
    deconv_2_kernel = graph.get_tensor_by_name('adeconv_2/kernel:0')

    deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                             2 * (code_size[3] - 1) + 5, 16])

    hidden_3_dec = tf.nn.conv3d_transpose(batch_n_d_2, deconv_2_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_2_bias

    hidden_3_dec = tf.nn.relu(hidden_3_dec)
    mu_3 = graph.get_tensor_by_name('abatch_n_d_2/moving_mean:0')
    var_3 = graph.get_tensor_by_name('abatch_n_d_2/moving_variance:0')
    beta_3 = graph.get_tensor_by_name('abatch_n_d_2/beta:0')
    gamma_3 = graph.get_tensor_by_name('abatch_n_d_2/gamma:0')
    batch_n_d_3 = gamma_3 * ((hidden_3_dec - mu_3) / tf.sqrt(var_3 + epsilon)) + beta_3

    # Fourth deconv layer
    code_size = tf.shape(batch_n_d_3)
    deconv_3_bias = graph.get_tensor_by_name('adeconv_3/bias:0')
    deconv_3_kernel = graph.get_tensor_by_name('adeconv_3/kernel:0')

    deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                             2 * (code_size[3] - 1) + 5, 24])

    hidden_4_dec = tf.nn.conv3d_transpose(batch_n_d_3, deconv_3_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_3_bias

    hidden_4_dec = tf.nn.relu(hidden_4_dec)
    mu_4 = graph.get_tensor_by_name('abatch_n_d_3/moving_mean:0')
    var_4 = graph.get_tensor_by_name('abatch_n_d_3/moving_variance:0')
    beta_4 = graph.get_tensor_by_name('abatch_n_d_3/beta:0')
    gamma_4 = graph.get_tensor_by_name('abatch_n_d_3/gamma:0')
    batch_n_d_4 = gamma_4 * ((hidden_4_dec - mu_4) / tf.sqrt(var_4 + epsilon)) + beta_4

    # Fifth deconv layer
    code_size = tf.shape(batch_n_d_4)
    deconv_4_bias = graph.get_tensor_by_name('adeconv_4/bias:0')
    deconv_4_kernel = graph.get_tensor_by_name('adeconv_4/kernel:0')

    deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                             2 * (code_size[3] - 1) + 5, 32])

    hidden_5_dec = tf.nn.conv3d_transpose(batch_n_d_4, deconv_4_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_4_bias

    hidden_5_dec = tf.nn.relu(hidden_5_dec)
    mu_5 = graph.get_tensor_by_name('abatch_n_d_4/moving_mean:0')
    var_5 = graph.get_tensor_by_name('abatch_n_d_4/moving_variance:0')
    beta_5 = graph.get_tensor_by_name('abatch_n_d_4/beta:0')
    gamma_5 = graph.get_tensor_by_name('abatch_n_d_4/gamma:0')
    batch_n_d_5 = gamma_5 * ((hidden_5_dec - mu_5) / tf.sqrt(var_5 + epsilon)) + beta_5

    # Sixth deconv layer
    code_size = tf.shape(batch_n_d_5)
    deconv_5_bias = graph.get_tensor_by_name('adeconv_5/bias:0')
    deconv_5_kernel = graph.get_tensor_by_name('adeconv_5/kernel:0')
    deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                             2 * (code_size[3] - 1) + 5, 1])

    hidden_6_dec = tf.nn.conv3d_transpose(batch_n_d_5, deconv_5_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_5_bias

    # We put -100 to the padding to be sure that the final prob is almost if not zero
    hidden_6_dec = pad_up_to(hidden_6_dec, [1, imageSize[0], imageSize[1], imageSize[2], 1],
                             constant_values=-100)

    return tf.nn.sigmoid(hidden_6_dec)


# add paddings when the size last layer does not match the input size, this happens for deconvolution layers
def pad_up_to(t, max_in_dims, constant_values):
    s = tf.shape(t)
    paddings = [[tf.cast((m - s[i]) / 2, tf.int32), (m - s[i]) - tf.cast((m - s[i]) / 2, tf.int32)]
                for (i, m) in enumerate(max_in_dims)]
    return tf.pad(t, tf.cast(tf.stack(paddings), tf.int32), 'CONSTANT', constant_values=constant_values)


# here is the net function. The input goes through the encoder, we sample from it and then it goes through the decoder
def run_net(lesion, imageSize):
    mu, sigma = get_encoder(lesion)
    sample_latent = tf.random.normal(mu.shape, 0, 1) * sigma + mu
    return get_decoder(sample_latent, imageSize)


def checkConditions(searchString, checkStructureOwnClass=True):

    global atlasDir_

    # The implementation here only works for the special scenario where
    #   (a) The structure of interest has its own class (mixture model) not shared with any other structure
    #       Not checked if checkStructureOwnClass=False
    #   (b) This class (mixture model) has only a single component
    #   (c) The structure of interest is not a mixture of two or more classes (mixture models)
    # Let's test for that here

    # Get class number
    classNumber = getClassNumber(searchString)

    # Get class fractions
    modelSpecifications = getModelSpecifications(atlasDir_)
    modelSpecifications = Specification(modelSpecifications)
    numberOfGaussiansPerClass = [ param.numberOfComponents for param in modelSpecifications.sharedGMMParameters ]

    classFractions, _ = gems.kvlGetMergingFractionsTable( modelSpecifications.names,
                                                          modelSpecifications.sharedGMMParameters )

    structureNumbers = np.flatnonzero(classFractions[classNumber, :] == 1)
    gaussianNumbers = [sum(numberOfGaussiansPerClass[0: classNumber])]
    if checkStructureOwnClass and structureNumbers.size is not 1:
        raise Exception('Structure of interest should correspond to exactly one class (mixture model) and vice versa')
    if len(gaussianNumbers) is not 1:
        raise Exception('Structure of interest should have a mixture model with only a single component')

    return structureNumbers[0], classNumber, gaussianNumbers[0]


def evaluateMinLogPriorOfGMMParametersWithLesCost(means, variances, mixtureWeights, numberOfGaussiansPerClass,
                                                  hyperMeans, hyperMeansNumberOfMeasurements,
                                                  hyperVariances, hyperVariancesNumberOfMeasurements,
                                                  hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements,
                                                  evaluateMinLogPriorPluginDictionary=None):

    global lesionGaussianNumber_, wmGaussianNumber_, numberOfPseudoSamplesMean_, numberOfPseudoSamplesVariance_, rho_

    # Make sure the hyperparameters are defined and valid
    numberOfContrasts = means.shape[1]
    hyperMeans, hyperMeansNumberOfMeasurements, \
           hyperVariances, hyperVariancesNumberOfMeasurements, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements = \
               getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts,
                                        hyperMeans, hyperMeansNumberOfMeasurements,
                                        hyperVariances, hyperVariancesNumberOfMeasurements,
                                        hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements )

    # We need to take into account the changes in the hyperprior of the lesion gaussian that is tied to the wm one
    hyperMeans[lesionGaussianNumber_] = means[wmGaussianNumber_]
    hyperMeansNumberOfMeasurements[lesionGaussianNumber_] = numberOfPseudoSamplesMean_
    hyperVariances[lesionGaussianNumber_] = rho_ * variances[wmGaussianNumber_]
    hyperVariancesNumberOfMeasurements[lesionGaussianNumber_] = numberOfPseudoSamplesVariance_

    # Now call the default evaluateMinLogPriorOfGMMParameters function with the updated parameters
    minLogGMMParametersPrior = evaluateMinLogPriorOfGMMParameters(means, variances, mixtureWeights,
                                                                  numberOfGaussiansPerClass,
                                                                  hyperMeans,
                                                                  hyperMeansNumberOfMeasurements,
                                                                  hyperVariances,
                                                                  hyperVariancesNumberOfMeasurements,
                                                                  hyperMixtureWeights,
                                                                  hyperMixtureWeightsNumberOfMeasurements)

    return minLogGMMParametersPrior


def fitGMMLesWMTied(data, gaussianPosteriors, numberOfGaussiansPerClass, useDiagonalCovarianceMatrices,
                    hyperMeans, hyperMeansNumberOfMeasurements, hyperVariances, hyperVariancesNumberOfMeasurements,
                    hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements, fitGMMPluginDictionary):

    global lesionGaussianNumber_, wmGaussianNumber_, numberOfPseudoSamplesMean_, numberOfPseudoSamplesVariance_,  rho_

    # First call the default fitGMMParameters first, then we "overwrite" the means and variances of lesion and WM
    means, variances, mixtureWeights = fitGMMParameters(data, gaussianPosteriors,
                                                        numberOfGaussiansPerClass,
                                                        useDiagonalCovarianceMatrices,
                                                        hyperMeans,
                                                        hyperMeansNumberOfMeasurements,
                                                        hyperVariances,
                                                        hyperVariancesNumberOfMeasurements,
                                                        hyperMixtureWeights,
                                                        hyperMixtureWeightsNumberOfMeasurements)

    # Make sure the hyperparameters are defined and valid
    numberOfContrasts = data.shape[1]
    hyperMeans, hyperMeansNumberOfMeasurements, \
           hyperVariances, hyperVariancesNumberOfMeasurements, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements = \
               getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts,
                                        hyperMeans, hyperMeansNumberOfMeasurements,
                                        hyperVariances, hyperVariancesNumberOfMeasurements,
                                        hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements )

    # Now update the means and variances of the wm and lesion class
    previousVariances = fitGMMPluginDictionary['variances']

    posterior_wm = gaussianPosteriors[:, wmGaussianNumber_].reshape(-1, 1)
    hyperMean_wm = np.expand_dims(hyperMeans[wmGaussianNumber_, :], 1)
    hyperMeanNumberOfMeasurements_wm = hyperMeansNumberOfMeasurements[wmGaussianNumber_]
    hyperVariance_wm = hyperVariances[wmGaussianNumber_, :, :]
    hyperVarianceNumberOfMeasurements_wm = hyperVariancesNumberOfMeasurements[wmGaussianNumber_]
    variance_wm_previous = previousVariances[wmGaussianNumber_]

    posterior_les = gaussianPosteriors[:, lesionGaussianNumber_].reshape(-1, 1)
    hyperMeanNumberOfMeasurements_les = numberOfPseudoSamplesMean_
    hyperVarianceNumberOfMeasurements_les = numberOfPseudoSamplesVariance_
    variance_les_previous = previousVariances[lesionGaussianNumber_]

    # Define some temporary variables
    soft_sum_wm = np.sum(posterior_wm)
    soft_sum_les = np.sum(posterior_les)
    wmLesTiedMean = ((hyperMeanNumberOfMeasurements_les * soft_sum_les) /
                     (soft_sum_les + hyperMeanNumberOfMeasurements_les)) * \
                    np.linalg.solve(variance_les_previous, variance_wm_previous)

    # Updates for the mean of wm and lesion gaussians
    mean_wm = np.linalg.solve((soft_sum_wm + hyperMeanNumberOfMeasurements_wm) * np.eye(numberOfContrasts) + wmLesTiedMean,
                              (data.T @ posterior_wm + hyperMean_wm * hyperMeanNumberOfMeasurements_wm +
                               wmLesTiedMean @ ((data.T @ posterior_les) / soft_sum_les)))

    mean_les = (data.T @ posterior_les + mean_wm * hyperMeanNumberOfMeasurements_les) /\
               (soft_sum_les + hyperMeanNumberOfMeasurements_les)

    # Updates for the variance of wm and lesion gaussians
    tmp = data - mean_wm.T

    wmLesTiedVariance = (2 + numberOfContrasts) * np.eye(numberOfContrasts) + hyperVarianceNumberOfMeasurements_les *\
                        (np.linalg.solve(variance_les_previous, rho_ * variance_wm_previous) - np.eye(numberOfContrasts))

    variance_wm = np.linalg.solve((soft_sum_wm + hyperVarianceNumberOfMeasurements_wm) *
                                  np.eye(numberOfContrasts) + wmLesTiedVariance,
                                  (tmp.T @ (tmp * posterior_wm) + hyperMeanNumberOfMeasurements_wm *
                                   ((mean_wm - hyperMean_wm) @ (mean_wm - hyperMean_wm).T) + hyperVariance_wm))

    tmp = data - mean_les.T

    variance_les = (tmp.T @ (tmp * posterior_les) + hyperMeanNumberOfMeasurements_les *
                    ((mean_les - mean_wm) @ (mean_les - mean_wm).T) +
                    hyperVarianceNumberOfMeasurements_les * rho_ * variance_wm) \
                   / (soft_sum_les + hyperVarianceNumberOfMeasurements_les)

    if useDiagonalCovarianceMatrices:
        variance_les = np.diag(np.diag(variance_les))
        variance_wm = np.diag(np.diag(variance_wm))

    means[wmGaussianNumber_, :] = mean_wm.T
    means[lesionGaussianNumber_, :] = mean_les.T
    variances[wmGaussianNumber_, :, :] = variance_wm
    variances[lesionGaussianNumber_, :, :] = variance_les

    return means, variances, mixtureWeights


def getMCMCPosteriors(data, priors, means, variances, mixtureWeights,
                      numberOfGaussiansPerClass, classFractions,
                      posteriorPluginDictionary):

    # Retrieve the global variables
    global atlasDir_, lesionStructureNumber_, lesionClassNumber, lesionGaussianNumber_, numberOfPseudoSamplesMean_,\
        numberOfPseudoSamplesVariance_, rho_, wmGaussianNumber_, intensityMaskingClassNumber_,\
        numberOfSamplingSteps_, numberOfBurnInSteps_, intensityMaskingPattern_, visualizer_

    # Unpack the extra variables we asked for
    useDiagonalCovarianceMatrices = posteriorPluginDictionary['modelSpecifications.useDiagonalCovarianceMatrices']
    mask = posteriorPluginDictionary['mask']
    voxelSpacing = posteriorPluginDictionary['voxelSpacing']
    transform = posteriorPluginDictionary['transform']

    #
    numberOfClasses = len(numberOfGaussiansPerClass)
    numberOfVoxels = data.shape[0]
    numberOfContrasts = data.shape[-1]
    imageSize = mask.shape

    # Create intensity-based lesion mask
    if intensityMaskingClassNumber_ is not None:
        # We have -1 mask below mean, +1 above, 0 nothing
        intensityMaskingGaussianNumbers = [sum(numberOfGaussiansPerClass[0: intensityMaskingClassNumber_])]
        intensityMaskingMean = np.sum(means[intensityMaskingGaussianNumbers] *
                                      mixtureWeights[intensityMaskingGaussianNumbers].reshape(-1, 1), axis=0)
        tmp = np.ones(numberOfVoxels, dtype=bool)
        for contrastNumber in range(numberOfContrasts):
            direction = intensityMaskingPattern_[contrastNumber]
            if direction == -1:
                tmp = np.logical_and(tmp, data[:, contrastNumber] < intensityMaskingMean[contrastNumber])
            elif direction == 1:
                tmp = np.logical_and(tmp, data[:, contrastNumber] > intensityMaskingMean[contrastNumber])
        intensityMask = np.zeros(imageSize, dtype=bool)
        intensityMask[mask] = tmp
    else:
        intensityMask = np.zeros(imageSize, dtype=bool)
        intensityMask[mask] = True

    visualizer_.show(image_list=[intensityMask.astype(float)], title="Intensity mask")

    # Initialize the structure likelihoods from the initial parameter values. Since only the parameters of a single structure
    # will be altered, only one column in the likelihoods will need to be updated during sampling
    likelihoods = getLikelihoods(data, means, variances, mixtureWeights, numberOfGaussiansPerClass, classFractions)

    # Initialize the sampler with a majority-vote lesion segmentation, masked with intensity-based mask
    posteriors = likelihoods * priors
    posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
    lesion = np.zeros(imageSize)
    lesion[mask] = (np.array(np.argmax(posteriors, 1), dtype=np.uint32) == lesionStructureNumber_)
    lesion *= intensityMask

    visualizer_.show(image_list=[lesion], title="Initial lesion segmentation")

    # Initalize the VAE tensorflow model and its various settings.
    # TODO: remove various warnings about deprecated function calls etc

    # Restore from checkpoint
    sess = tf.Session()
    VAEModelPath = os.path.join(atlasDir_, 'VAE')
    VAEModelFileName = os.path.join(VAEModelPath, 'model.ckpt.meta')
    saver = tf.train.import_meta_graph(VAEModelFileName)
    saver.restore(sess, tf.train.latest_checkpoint(VAEModelPath))
    print('VAE lesion model loaded')

    # Get some info about the VAE
    VAEInfo = np.load(os.path.join(atlasDir_, 'VAE/VAEInfo.npz'))
    net_shape = VAEInfo['net_shape']
    # transformation matrix mapping train VAE coordinates to template coordinates
    trainToTemplateMat = VAEInfo['trainToTemplateMat']

    # Combination of transformation matrixes in order to obtain a subject to VAE train space transformation
    # First from subject space to template space, then from template space to VAE train space
    # When combining transformations the order of the transformations is from right to left.
    trainToSubjectMat = transform.as_numpy_array @ trainToTemplateMat
    subjectToTrainMat = np.linalg.inv(trainToSubjectMat)

    # Create tf placeholder
    lesionPlaceholder = tf.placeholder(tf.float32, [1, net_shape[0], net_shape[1], net_shape[2], 1])
    net = run_net(lesionPlaceholder, net_shape)

    # Do the actual sampling of lesion, latent variables of the VAE model, and mean/variance of the lesion intensity model.
    hyperMean = means[wmGaussianNumber_]
    hyperMeanNumberOfMeasurements = numberOfPseudoSamplesMean_
    hyperVariance = rho_ * variances[wmGaussianNumber_]
    hyperVarianceNumberOfMeasurements = numberOfPseudoSamplesVariance_
    averagePosteriors = np.zeros_like(likelihoods)
    visualizer_.start_movie(window_id="Lesion prior using VAE only", title="Lesion prior using VAE only -- the movie")
    visualizer_.start_movie(window_id="Lesion sample", title="Lesion sample -- the movie")
    for sweepNumber in range(numberOfBurnInSteps_ + numberOfSamplingSteps_):

        # Sample from the VAE latent variables, conditioned on the current lesion segmentation.
        # Implementation-wise we don't store the latent variables, but rather the factorized
        # prior in the visible units (voxels) that they encode.

        # We first go from subject space to train space of the VAE
        # Since we are using scipy affine transform that takes an INVERSE transformation
        # we pass to the function the inverse of subjectToTrainMat, so trainToSubjectMat
        lesionTrainSpace = affine_transform(lesion, trainToSubjectMat, output_shape=net_shape, order=1)

        # We go through the VAE to get the factorized prior
        # tensorflow wants a 5 dimensional input where the first dimension is batch number
        # and last dimension is number of channels, both set to 1 in our case.
        expandedLesion = np.expand_dims(np.expand_dims(lesionTrainSpace, 0), 4)
        lesionVAETrainSpace = np.squeeze(np.squeeze(sess.run(net, {lesionPlaceholder: expandedLesion}), 4), 0)

        # We then go back to subject space from train space
        # Also here, since we are using scipy affine transform that takes an INVERSE transformation
        # we pass to the function the inverse of trainToSubjectMat, so subjectToTrainMat
        lesionPriorVAE = (affine_transform(lesionVAETrainSpace, subjectToTrainMat, output_shape=imageSize, order=1)
                          * intensityMask)[mask]

        if hasattr(visualizer_, 'show_flag'):
            tmp = np.zeros(imageSize)
            tmp[mask] = lesionPriorVAE
            visualizer_.show(probabilities=tmp, title="Lesion prior using VAE only",
                             window_id="Lesion prior using VAE only")

        # Sample from the mean and variance, conditioned on the data and the lesion segmentation
        bestMean, bestVariance, _ = fitGMMParameters(data,
                                                     lesion[mask].reshape(-1, 1),
                                                     [1], useDiagonalCovarianceMatrices,
                                                     hyperMeans=np.array([hyperMean]),
                                                     hyperMeansNumberOfMeasurements=np.array(
                                                         [hyperMeanNumberOfMeasurements]),
                                                     hyperVariances=np.array([hyperVariance]),
                                                     hyperVariancesNumberOfMeasurements=np.array(
                                                         [hyperVarianceNumberOfMeasurements])
                                                     )
        bestMean, bestVariance = bestMean[0], bestVariance[0]
        N = lesion[mask].sum()

        # Murphy, page 134 with v0 = hyperVarianceNumberOfMeasurements - numberOfContrasts - 2
        variance = invwishart.rvs(N + hyperVarianceNumberOfMeasurements - numberOfContrasts - 2,
                                  bestVariance * (hyperVarianceNumberOfMeasurements + N))

        # If numberOfContrast is 1 force variance to be a (1,1) array
        if numberOfContrasts == 1:
            variance = np.atleast_2d(variance)

        if useDiagonalCovarianceMatrices:
            variance = np.diag(np.diag(variance))

        mean = np.random.multivariate_normal(bestMean, variance / (hyperMeanNumberOfMeasurements + N)).reshape(-1, 1)

        # Sample from the lesion segmentation, conditioned on the data and the VAE latent variables
        # (Implementation-wise the latter is encoded in the VAE prior). At the same time we also
        # compute the full posterior of each structure, which is at the end the thing we're averaging
        # over (i.e., the reason why we're sampling)
        effectivePriors = np.array(priors / 65535, dtype=np.float32)
        effectivePriors[:, lesionStructureNumber_] *= lesionPriorVAE
        if hasattr(visualizer_, 'show_flag'):
            tmp = np.zeros(imageSize)
            tmp[mask] = effectivePriors[:, lesionStructureNumber_]
            visualizer_.show(probabilities=tmp, title="Lesion prior using VAE and atlas together",
                             window_id="Lesion prior using VAE and atlas together")

        # Generative model where the atlas generates *candidate* lesions, and the VAE prior is sampled
        # from *only within the candidates*.
        numberOfStructures = priors.shape[-1]
        otherStructureNumbers = [i for i in range(numberOfStructures) if i != lesionStructureNumber_]
        otherStructurePriors = effectivePriors[:, otherStructureNumbers]
        otherStructurePriors /= (np.sum(otherStructurePriors, axis=1).reshape(-1, 1) + eps)
        effectivePriors[:, otherStructureNumbers] = otherStructurePriors * \
                                                    (1 - effectivePriors[:, lesionStructureNumber_].reshape(-1, 1))
        likelihoods[:, lesionStructureNumber_] = getGaussianLikelihoods(data, mean, variance)
        posteriors = effectivePriors * likelihoods
        posteriors /= np.expand_dims(np.sum(posteriors, axis=1) + eps, 1)
        sample = np.random.rand(numberOfVoxels) <= posteriors[:, lesionStructureNumber_]
        lesion = np.zeros(imageSize)
        lesion[mask] = sample

        visualizer_.show(image_list=[lesion], title="Lesion sample", window_id="Lesion sample")

        # Collect data after burn in steps
        if sweepNumber >= numberOfBurnInSteps_:
            print('Sample ' + str(sweepNumber + 1 - numberOfBurnInSteps_) + ' times')
            averagePosteriors += posteriors / numberOfSamplingSteps_
        else:
            print('Burn-in ' + str(sweepNumber + 1) + ' times')

    #
    visualizer_.show_movie(window_id="Lesion prior using VAE only")
    visualizer_.show_movie(window_id="Lesion sample")

    # Return
    return averagePosteriors


def samsegmentLesion(imageFileNames, atlasDir, savePath,
                     transformedTemplateFileName=None,
                     userModelSpecifications={}, userOptimizationOptions={},
                     lesionSearchString='Lesion',
                     numberOfSamplingSteps=50, numberOfBurnInSteps=50,
                     visualizer=None, saveHistory=False, saveMesh=False,
                     saveWarp=False, targetIntensity=None, targetSearchStrings=None,
                     intensityMaskingPattern=None, intensityMaskingSearchString=None,
                     lesionThreshold=None,
                     savePosteriors=False,
                     numberOfPseudoSamplesMean=500,
                     numberOfPseudoSamplesVariance=500,
                     rho=50
                     ):
    # We can't pass on user-specified options to plugins, so let's do it using global variables
    # (using underscore notation to stress they're global variables )
    global numberOfPseudoSamplesMean_, numberOfPseudoSamplesVariance_, rho_, atlasDir_, lesionStructureNumber_,\
        lesionClassNumber_, lesionGaussianNumber_, wmGaussianNumber_, intensityMaskingClassNumber_,\
        numberOfSamplingSteps_, numberOfBurnInSteps_, intensityMaskingPattern_, visualizer_

    atlasDir_ = atlasDir

    # Check conditions on white matter and lesion gaussian/structure and
    # get their structure numbers, class number as well as the gaussian number
    wmSearchString = 'White'
    lesionStructureNumber_, lesionClassNumber_, lesionGaussianNumber_ = checkConditions(lesionSearchString)
    _, _, wmGaussianNumber_ = checkConditions(wmSearchString, checkStructureOwnClass=False)

    numberOfPseudoSamplesMean_ = numberOfPseudoSamplesMean
    numberOfPseudoSamplesVariance_ = numberOfPseudoSamplesVariance
    rho_ = rho
    intensityMaskingClassNumber_ = getClassNumber(intensityMaskingSearchString)
    numberOfSamplingSteps_ = numberOfSamplingSteps
    numberOfBurnInSteps_ = numberOfBurnInSteps
    intensityMaskingPattern_ = intensityMaskingPattern
    if visualizer is None:
        visualizer = initVisualizer(False, False)
    visualizer_ = visualizer

    # Now call samsegment with plugins
    samsegment(imageFileNames, atlasDir, savePath,
               userModelSpecifications=userModelSpecifications,
               userOptimizationOptions=userOptimizationOptions,
               transformedTemplateFileName=transformedTemplateFileName,
               visualizer=visualizer,
               targetIntensity=targetIntensity,
               targetSearchStrings=targetSearchStrings,
               evaluateMinLogPriorOfGMMParametersPlugin=evaluateMinLogPriorOfGMMParametersWithLesCost,
               evaluateMinLogPriorOfGMMParametersPluginVariables=None,
               fitGMMParametersPlugin=fitGMMLesWMTied,
               fitGMMParametersPluginVariables=['variances'],
               posteriorPlugin=getMCMCPosteriors,
               posteriorPluginVariables=['modelSpecifications.useDiagonalCovarianceMatrices',
                                         'mask', 'voxelSpacing', 'transform'],
               threshold=lesionThreshold, thresholdSearchString=lesionSearchString,
               savePosteriors=savePosteriors, saveHistory=saveHistory,
               saveMesh=saveMesh, saveWarp=saveWarp
               )
