from freesurfer.samseg.samsegment import samsegment, initializeGMMParameters, getFullHyperparameters, \
                                         getLikelihoods, fitGMMParameters, getGaussianLikelihoods
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow.contrib.distributions as tfd
from scipy.stats import invwishart
from scipy.ndimage.interpolation import zoom
from freesurfer.gems import kvlReadCompressionLookupTable, kvlReadSharedGMMParameters
import os

#eps = np.spacing(1)
eps = np.finfo( float ).eps


  


def getClassNumber( structureSearchString ):
    #
    if structureSearchString is None:
        return None
    
    #
    global atlasDir_
  
    sharedGMMParameters = kvlReadSharedGMMParameters( os.path.join( atlasDir_, 'sharedGMMParameters.txt') )

    for classNumber, mergeOption in enumerate( sharedGMMParameters ):
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

    return tfd.MultivariateNormalDiag(mu, sigma)


# decoder network
def get_decoder(code, imageSize):
    epsilon = 1e-3
    
    graph = tf.get_default_graph()
    
    code_size = tf.shape(code)

    # First deconv layer
    deconv_0_bias = graph.get_tensor_by_name('adeconv_0/bias:0')
    deconv_0_kernel = graph.get_tensor_by_name('adeconv_0/kernel:0')
    # Deconvolution shape for VALID = stride * (input - 1) + kernel size
    deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3, code_size[3] - 1 + 3, 16])

    hidden_0_dec = tf.nn.conv3d_transpose(code, deconv_0_kernel, output_shape=deconv_shape,
                                          strides=[1, 1, 1, 1, 1],
                                          padding='VALID', ) + deconv_0_bias

    mu_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_mean:0')
    var_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_variance:0')
    beta_1 = graph.get_tensor_by_name('abatch_n_d_0/beta:0')
    gamma_1 = graph.get_tensor_by_name('abatch_n_d_0/gamma:0')
    batch_n_d_1 = gamma_1 * ((hidden_0_dec - mu_1) / tf.sqrt(var_1 + epsilon)) + beta_1

    # Second deconv layer
    code_size = tf.shape(batch_n_d_1)
    deconv_1_bias = graph.get_tensor_by_name('adeconv_1/bias:0')
    deconv_1_kernel = graph.get_tensor_by_name('adeconv_1/kernel:0')
    deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3, code_size[3] - 1 + 3, 16])

    hidden_2_dec = tf.nn.conv3d_transpose(batch_n_d_1, deconv_1_kernel, output_shape=deconv_shape,
                                          strides=[1, 1, 1, 1, 1],
                                          padding='VALID') + deconv_1_bias

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

    deconv_shape = tf.stack(
        [code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5, 2 * (code_size[3] - 1) + 5, 16])

    hidden_3_dec = tf.nn.conv3d_transpose(batch_n_d_2, deconv_2_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1],
                                          padding='VALID') + deconv_2_bias

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

    deconv_shape = tf.stack(
        [code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5, 2 * (code_size[3] - 1) + 5, 24])

    hidden_4_dec = tf.nn.conv3d_transpose(batch_n_d_3, deconv_3_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1],
                                          padding='VALID') + deconv_3_bias

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

    deconv_shape = tf.stack(
        [code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5, 2 * (code_size[3] - 1) + 5, 32])

    hidden_5_dec = tf.nn.conv3d_transpose(batch_n_d_4, deconv_4_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1],
                                          padding='VALID') + deconv_4_bias

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
    deconv_shape = tf.stack(
        [code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5, 2 * (code_size[3] - 1) + 5, 1])

    hidden_6_dec = tf.nn.conv3d_transpose(batch_n_d_5, deconv_5_kernel, output_shape=deconv_shape,
                                          strides=[1, 2, 2, 2, 1],
                                          padding='VALID') + deconv_5_bias

    # We put -100 to the padding to be sure that the final prob is almost if not zero
    hidden_6_dec = pad_up_to(hidden_6_dec, [1, imageSize[0], imageSize[1], imageSize[2], 1], constant_values=-100)

    return tf.nn.sigmoid(hidden_6_dec)


# add paddings when the size last layer does not match the input size, this happens for deconvolution layers
def pad_up_to(t, max_in_dims, constant_values):
    s = tf.shape(t)
    paddings = [[tf.cast((m - s[i]) / 2, tf.int32), (m - s[i]) - tf.cast((m - s[i]) / 2, tf.int32)] for (i, m) in
                enumerate(max_in_dims)]
    return tf.pad(t, tf.cast(tf.stack(paddings), tf.int32), 'CONSTANT', constant_values=constant_values)


# here is the net function. The input goes through the encoder, we sample from it and then it goes through the decoder
def run_net(lesion, imageSize):
    code = get_encoder(lesion).sample()
    return get_decoder(code, imageSize)


# pad function for scaling to image size from net size
def scale(data, imageSize, net_size):
    diff_1 = ((net_size - imageSize) / 2).astype(int)
    diff_2 = ((net_size - imageSize) - diff_1).astype(int)
    if diff_2[0] < 0:
        diff_1[0] = 0
        diff_2[0] = 0
    if diff_2[1] < 0:
        diff_1[1] = 0
        diff_2[1] = 0
    if diff_2[2] < 0:
        diff_1[2] = 0
        diff_2[2] = 0
    if len(data.shape) == 3:
        paddings = [[diff_1[0], diff_2[0]], [diff_1[1], diff_2[1]], [diff_1[2], diff_2[2]]]
    else:
        paddings = [[0, 0], [diff_1[0], diff_2[0]], [diff_1[1], diff_2[1]], [diff_1[2], diff_2[2]], [0, 0]]

    return np.pad(data, paddings, mode='constant')


# pad function for cropping to image size from net size
def crop(data, imageSize, net_size):
    diff_1 = ((net_size - imageSize) / 2).astype(int)
    diff_2 = ((net_size - imageSize) - diff_1).astype(int)
    if diff_2[0] < 0:
        diff_1[0] = 0
        diff_2[0] = 0
    if diff_2[1] < 0:
        diff_1[1] = 0
        diff_2[1] = 0
    if diff_2[2] < 0:
        diff_1[2] = 0
        diff_2[2] = 0
    if len(data.shape) == 3:
        return data[diff_1[0]:(imageSize[0] + diff_1[0]), diff_1[1]:(imageSize[1] + diff_1[1]),
               diff_1[2]:(imageSize[2] + diff_1[2])]
    data = data[0, diff_1[0]:(imageSize[0] + diff_1[0]), diff_1[1]:(imageSize[1] + diff_1[1]),
           diff_1[2]:(imageSize[2] + diff_1[2])]
    return data



def defineHyperparameters( data, classPriors, numberOfGaussiansPerClass, voxelSpacing ):
    #
    global hyperMeans_, hyperMeansNumberOfMeasurements_, hyperVariances_, hyperVariancesNumberOfMeasurements_, \
           lesionClassNumber_
    
    #
    initialMeans, _, _ = initializeGMMParameters( data, classPriors, numberOfGaussiansPerClass )
    lesionGaussianNumbers = [ sum( numberOfGaussiansPerClass[ 0:lesionClassNumber_ ] ) ]
    lesionHyperMeans = initialMeans[ lesionGaussianNumbers ]
    if False:
        lesionHyperMeansNumberOfMeasurements = 100
        lesionHyperVariancesNumberOfMeasurements = 500
    else:
        lesionHyperMeansNumberOfMeasurements = 100 / np.prod( voxelSpacing )
        lesionHyperVariancesNumberOfMeasurements = 500 / np.prod( voxelSpacing )
        dataVariance = np.var( data, axis=0 )
        if False:
            # Woah divided by lesionHyperVariancesNumberOfMeasurements -- why??
            lesionHyperVariances = 200 * np.diag( dataVariance ) / lesionHyperVariancesNumberOfMeasurements
        else:
            lesionHyperVariances = 0.4 * np.diag( dataVariance )

    numberOfContrasts = data.shape[ -1 ]
    hyperMeans_, hyperMeansNumberOfMeasurements_, \
           hyperVariances_, hyperVariancesNumberOfMeasurements_, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements = \
             getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts )

    hyperMeans_[ lesionGaussianNumbers ] = lesionHyperMeans
    hyperMeansNumberOfMeasurements_[ lesionGaussianNumbers ] = lesionHyperMeansNumberOfMeasurements
    hyperVariances_[ lesionGaussianNumbers, :, : ] = lesionHyperVariances
    hyperVariancesNumberOfMeasurements_[ lesionGaussianNumbers ] = lesionHyperVariancesNumberOfMeasurements

    return hyperMeans_, hyperMeansNumberOfMeasurements_, \
           hyperVariances_, hyperVariancesNumberOfMeasurements_, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements



def getMCMCPosteriors( data, priors, means, variances, mixtureWeights, 
                         numberOfGaussiansPerClass, classFractions,
                         posteriorPluginDictionary ):

    # Retrieve the global variables
    global hyperMeans_, hyperMeansNumberOfMeasurements_, hyperVariances_, hyperVariancesNumberOfMeasurements_, \
           atlasDir_, lesionClassNumber_, intensityMaskingClassNumber_, numberOfSamplingSteps_, numberOfBurnInSteps_, \
           intensityMaskingPattern_, visualizer_  


    # Unpack the extra variables we asked for
    useDiagonalCovarianceMatrices = posteriorPluginDictionary[ 'modelSpecifications.useDiagonalCovarianceMatrices' ]
    mask = posteriorPluginDictionary[ 'mask' ]
    voxelSpacing = posteriorPluginDictionary[ 'voxelSpacing' ]
    transform = posteriorPluginDictionary[ 'transform' ]

    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfVoxels = data.shape[ 0 ]
    numberOfContrasts = data.shape[ -1 ]
    imageSize = mask.shape
    lesionGaussianNumbers = [ sum( numberOfGaussiansPerClass[ 0 : lesionClassNumber_ ] ) ]


    # The implementation here only works for the special scenario where 
    #   (a) The structure you're sampling from has its own class (mixture model) not shared with any other structure
    #   (b) This class (mixture model) has only a single component
    #   (c) The structure you're sampling from is not a mixture of two or more classes (mixture models)
    # Let's test for that here
    lesionStructureNumbers = np.flatnonzero( classFractions[ lesionClassNumber_, : ] == 1 )
    if lesionStructureNumbers.size is not 1:
        raise Exception( 'Structure of interest should correspond to exactly one class (mixture model) and vice versa' )
    if len( lesionGaussianNumbers ) is not 1:
        raise Exception( 'Structure of interest should have a mixture model with only a single component' )
    lesionStructureNumber = lesionStructureNumbers[ 0 ]
    lesionGaussianNumber = lesionGaussianNumbers[ 0 ]



    # Create intensity-based lesion mask
    if intensityMaskingClassNumber_ is not None:
        # We have -1 mask below mean, +1 above, 0 nothing
        intensityMaskingGaussianNumbers = [ sum( numberOfGaussiansPerClass[ 0 : intensityMaskingClassNumber_ ] ) ]
        intensityMaskingMean = np.sum( means[ intensityMaskingGaussianNumbers ] * 
                                      mixtureWeights[ intensityMaskingGaussianNumbers ].reshape( -1, 1 ), axis=0 )
        tmp = np.ones( numberOfVoxels, dtype=bool )
        for contrastNumber in range( numberOfContrasts ):
            direction = intensityMaskingPattern_[ contrastNumber ]
            if direction == '-1':
                tmp = np.logical_and( tmp, data[:, contrastNumber ] < intensityMaskingMean[ contrastNumber ] )
            elif direction == '1':
                tmp = np.logical_and( tmp, data[:, contrastNumber ] > intensityMaskingMean[ contrastNumber ] )
        intensityMask = np.zeros( imageSize, dtype=bool ); intensityMask[ mask ] = tmp
    else:
        intensityMask = np.zeros( imageSize, dtype=bool ); intensityMask[ mask ] = True

    visualizer_.show( image_list=[ intensityMask.astype( float ) ], title="Intensity mask" )


    # Initialize the structure likelihoods from the initial parameter values. Since only the parameters of a single structure
    # will be altered, only one column in the likelihoods will need to be updated during sampling
    likelihoods = getLikelihoods( data, means, variances, mixtureWeights, numberOfGaussiansPerClass, classFractions )


    # Initialize the sampler with a majority-vote lesion segmentation, masked with intensity-based mask
    posteriors = likelihoods * priors
    posteriors /= np.expand_dims( np.sum( posteriors, axis=1 ) + eps, 1 )
    lesion = np.zeros( imageSize )
    lesion[ mask ] = ( np.array( np.argmax( posteriors, 1 ), dtype=np.uint32 ) == lesionStructureNumber )
    lesion *= intensityMask

    visualizer_.show( image_list=[ lesion ], title="Initial lesion segmentation" )


    if False:
        maskIndices = mask 
        atlas = priors
        lesion_idx = lesionClassNumber_
        dataMask = data
        lesionsInit = np.where( lesion )
        gm_idx = intensityMaskingGaussianNumbers[ 0 ] 
        hyperMean = hyperMeans_[ lesionGaussianNumbers[ 0 ], : ] 
        hyperMeanNumberOfMeasurements = hyperMeansNumberOfMeasurements_[ lesionGaussianNumbers[ 0 ] ]
        hyperVariance = hyperVariances_[ lesionGaussianNumbers, :, : ]
        hyperVarianceNumberOfMeasurements = hyperVariancesNumberOfMeasurements_[ lesionGaussianNumbers[ 0 ] ]
        gaussianNumber_lesion = lesionGaussianNumbers[ 0 ] 
        fractionsTable = classFractions

        # Create intensit-based lesion mask
        # We have -1 mask below mean, +1 above, 0 nothing
        gm_means = means[ gm_idx, : ]
        intensityMask = np.ones( numberOfVoxels, dtype=bool )
        for contrastNumber in range( numberOfContrasts ):
            direction = intensityMaskingPattern_[ contrastNumber ]
            if direction == '-1':
                intensityMask = np.logical_and( intensityMask, data[:, contrastNumber ] < gm_means[ contrastNumber ] )
            elif direction == '1':
                intensityMask = np.logical_and( intensityMask, data[:, contrastNumber ] > gm_means[ contrastNumber ] )
        tmp = np.zeros( imageSize, dtype=bool ); tmp[ mask ] = intensityMask; intensityMask = tmp

        visualizer_.show( image_list=[ intensityMask.astype( float ) ]  )

        # Initialize lesion variable as intensity-maksed initial lesion segmentation
        lesion = lesion * intensityMask




    # Initalize the VAE tensorflow model and its various settings.
    # TODO: remove various warnings about deprecated function calls etc

    # Restore from checkpoint
    sess = tf.Session()
    VAEModelPath = os.path.join( atlasDir_, 'VAE' )
    VAEModelFileName = os.path.join( VAEModelPath, 'model.ckpt.meta' )
    saver = tf.train.import_meta_graph( VAEModelFileName )
    saver.restore( sess, tf.train.latest_checkpoint( VAEModelPath ) )
    print( 'VAE lesion model loaded' )

    if True:
        #
        # Compute some parameters to turn a 3D lesion image ("lesion") into the correct format for tensorflow
        #
        # TODO: the variables croppedLesion, zoomedLesion, expandedLesion and scaledLesion are not actually
        #       used in the remainder, so this code should be cleaned up. If it's difficult to define/compute
        #       the actual tensorflow-related parameters without their help, consider initializing them to
        #       None, and compute them inside the main loop by checking if they are still None. This will
        #       allow these things to be computed once, without any code replication (which is hard to 
        #       maintain)
        #
        # TODO: make this more Python-like. E.g., dilations = list( 1/voxelSpacing ) or perhaps even
        #       dilations = 1/voxelSpacing instead of writing
        #       dilations = [1 / voxelSpacing[0], 1 / voxelSpacing[1], 1 / voxelSpacing[2]]
        #
        
        # Center (crop??) image
        # TODO: why is this needed? Perhaps document the reason, or otherwise remove this entire cropping step
        coords = np.argwhere( mask )
        x0, y0, z0 = coords.min(axis=0)
        x1, y1, z1 = coords.max(axis=0) + 1  # slices are exclusive at the top
        croppedIntensityMask = intensityMask[ x0:x1, y0:y1, z0:z1 ]

        # Get voxel resolution
        dilations = [1 / voxelSpacing[0], 1 / voxelSpacing[1], 1 / voxelSpacing[2]]

        # Paddings for go back to original image size
        paddings = [ [ x0, imageSize[0] - x1 ], [ y0, imageSize[1] - y1 ], [ z0, imageSize[2] - z1 ] ]

        # Create lesion input
        croppedLesion = lesion[ x0:x1, y0:y1, z0:z1 ]
        croppedLesion = croppedLesion.astype( float )

        # Zoom it 1x1x1 mm resolution
        # TODO: unless rotations of e.g., 180 degrees of the image orientation is not any problem 
        #       at all for the model, you should resample the entire 3D lesion image to the 1mm 
        #       isotropic template space the VAE was apparently trained on. For this purpose you
        #       should use the transform matrix in the samseg code, which encodes how an image grid
        #       location in a MNI-like space maps into the image grid of the image being segmented.
        zoomedLesion = zoom( croppedLesion, voxelSpacing, order=1 )

        # If input image is smaller than net_shape we need to pad
        # Size of the training image for the VAE
        # TODO: this should not be hard-coded but rather read from the VAE directory 
        net_size = np.array( [ 197, 233, 189 ] )
        net_shape = np.array([ 197, 233, 189 ] )

        # TODO: this needs to be clean up and coded much more concisely. Use existing Python functionality where possible
        if net_shape[0] > zoomedLesion.shape[0] or net_shape[1] > zoomedLesion.shape[1] or net_shape[2] > zoomedLesion.shape[2]:
            pad = True
            if net_shape[0] < zoomedLesion.shape[0]:
                net_shape[0] = zoomedLesion.shape[0]
            if net_shape[1] < zoomedLesion.shape[1]:
                net_shape[1] = zoomedLesion.shape[1]
            if net_shape[2] < zoomedLesion.shape[2]:
                net_shape[2] = zoomedLesion.shape[2]
        else:
            net_shape = np.array( imageSize )
            pad = False
        #imageSize = np.array([imageSize[0], imageSize[1], imageSize[2]])

        # Create tf placeholder
        lesionPlaceholder = tf.placeholder(tf.float32, [1, net_shape[0], net_shape[1], net_shape[2], 1])
        net = run_net(lesionPlaceholder, net_shape)
        expandedLesion = np.expand_dims(np.expand_dims(zoomedLesion, 0), 4)

        if pad:
            scaledLesion = scale(expandedLesion, [expandedLesion.shape[1], expandedLesion.shape[2], expandedLesion.shape[3]], net_size)
        else:
            scaledLesion = expandedLesion



    # Do the actual sampling of lesion, latent variables of the VAE model, and mean/variance of the lesion intensity model.
    hyperMean = hyperMeans_[ lesionGaussianNumber ]
    hyperMeanNumberOfMeasurements = hyperMeansNumberOfMeasurements_[ lesionGaussianNumber ]
    hyperVariance = hyperVariances_[ lesionGaussianNumber ]
    hyperVarianceNumberOfMeasurements = hyperVariancesNumberOfMeasurements_[ lesionGaussianNumber ]
    averagePosteriors = np.zeros_like( likelihoods )
    visualizer_.start_movie( window_id="Lesion prior using VAE only", title="Lesion prior using VAE only -- the movie" )
    visualizer_.start_movie( window_id="Lesion sample", title="Lesion sample -- the movie" )
    for sweepNumber in range( numberOfBurnInSteps_ + numberOfSamplingSteps_ ):
        if False:
            # Sample mean and variance for lesion gaussian
            # First we sample from the means given the variances
            # then we sample from the variances given the means
            # we then recompute the likelihood for the lesion gaussian
            for t in range(numberOfClasses):
                numberOfComponents = numberOfGaussiansPerClass[t]
                for componentNumber in range(numberOfComponents):
                    gaussianNumber = int(np.sum(numberOfGaussiansPerClass[: t]) + componentNumber)
                    if gaussianNumber == gaussianNumber_lesion:
                        posterior = posteriors[:, lesion_idx]
                        posterior = posterior.reshape(-1, 1)

                        mean = (dataMask.T @ posterior + hyperMean.T * hyperMeanNumberOfMeasurements) \
                                / (np.sum(posterior) + hyperMeanNumberOfMeasurements)
                        means[gaussianNumber, :] = np.random.multivariate_normal(mean.ravel(),
                                                                                  variances[gaussianNumber, :, :] / (np.sum(
                                                                                      posterior) + hyperMeanNumberOfMeasurements))
                        tmp = dataMask - mean.T
                        S = tmp.T @ (tmp * posterior) + \
                            ((hyperMeanNumberOfMeasurements * np.sum(posterior)) / (np.sum(posterior) + hyperMeanNumberOfMeasurements)
                              * ((mean - hyperMean) @ (mean - hyperMean).T)) + \
                                hyperVariance * hyperVarianceNumberOfMeasurements
                        variances[gaussianNumber, :, :] = invwishart.rvs(
                            np.sum(posterior) + 1 + hyperVarianceNumberOfMeasurements, S)
                        if useDiagonalCovarianceMatrices:
                            # Force diagonal covariance matrices
                            variances[gaussianNumber, :, :] = np.diag(np.diag(variances[gaussianNumber, :, :]))
                        mean = np.expand_dims(means[gaussianNumber, :], 1)
                        variance = variances[gaussianNumber, :, :]

                        likelihoods[:, lesion_idx] = getGaussianLikelihoods( dataMask, mean, variance )
                        
                        
            # Run network and get prior
            if pad:
                data = crop(np.squeeze(np.squeeze(sess.run(net, {lesionPlaceholder: scaledLesion}), 4), 0),
                            [expandedLesion.shape[1], expandedLesion.shape[2], expandedLesion.shape[3]], net_size)
                data = np.clip(zoom(data, dilations), 0, 1)
                # If we have some rounding problems from zoom, pad with 0 for the rounding errors.
                if data.shape < croppedIntensityMask.shape:
                    data = scale(data, np.array([data.shape[0], data.shape[1], data.shape[2]]),
                                  np.array([croppedIntensityMask.shape[0], croppedIntensityMask.shape[1], croppedIntensityMask.shape[2]]))
                if data.shape > croppedIntensityMask.shape:
                    data = crop(data, np.array([croppedIntensityMask.shape[0], croppedIntensityMask.shape[1], croppedIntensityMask.shape[2]]),
                                np.array([data.shape[0], data.shape[1], data.shape[2]]))
            else:
                data = np.squeeze(np.squeeze(sess.run(net, {lesionPlaceholder: scaledLesion}), 4), 0)
                data = np.clip(zoom(data, dilations), 0, 1)

            # Assign lesion prior
            priors = np.array(atlas / 65535, dtype=np.float32)
            lesionPriorVAE = np.reshape(np.pad((data * croppedIntensityMask), paddings, mode='constant')[maskIndices == 1], -1).astype(
                np.float32)
            priors[:, lesion_idx] = lesionPriorVAE * priors[:, lesion_idx]
            normalizer = np.sum(priors, axis=1) + eps
            priors = priors / np.expand_dims(normalizer, 1)

            posteriors = priors * likelihoods

            # Normalize posteriors
            normalizer = np.sum(posteriors, axis=1) + eps
            posteriors = posteriors / np.expand_dims(normalizer, 1)

            # Multinomial sampling for posteriors in order to get new lesion sample
            print('Multinomial sampling')
            diff = np.cumsum(posteriors, axis=1) - np.random.random([posteriors.shape[0], 1])
            diff[diff < 0] = 1.1
            diff = np.argmin(diff, axis=1)
            lesion = np.zeros(imageSize[0] * imageSize[1] * imageSize[2])
            lesion[np.reshape(maskIndices == 1, -1)] = (diff == lesion_idx)
            lesion = np.reshape(lesion, imageSize)
            
            # Prepare data to feed to the network for next loop cycle
            croppedLesion = lesion[x0:x1, y0:y1, z0:z1]
            zoomedLesion = zoom(croppedLesion, voxelSpacing, order=1)
            expandedLesion = np.expand_dims(np.expand_dims(zoomedLesion, 0), 4)
            if pad:
                scaledLesion = scale(expandedLesion, [expandedLesion.shape[1], expandedLesion.shape[2], expandedLesion.shape[3]], net_size)
            else:
                scaledLesion = expandedLesion
                            
                            
        else:
            
            # Sample from the VAE latent variables, conditioned on the current lesion segmentation. 
            # Implementation-wise we don't store the latent variables, but rather the factorized
            # prior in the visible units (voxels) that they encode.
            # TODO: clean up the code, using more concise Python functionality and giving informative names 
            # to variables
            croppedLesion = lesion[ x0:x1, y0:y1, z0:z1 ]
            zoomedLesion = zoom(croppedLesion, voxelSpacing, order=1)
            expandedLesion = np.expand_dims(np.expand_dims(zoomedLesion, 0), 4)
            if pad:
                scaledLesion = scale(expandedLesion, [expandedLesion.shape[1], expandedLesion.shape[2], expandedLesion.shape[3]], net_size)
                zzz = crop(np.squeeze(np.squeeze(sess.run(net, {lesionPlaceholder: scaledLesion}), 4), 0),
                            [expandedLesion.shape[1], expandedLesion.shape[2], expandedLesion.shape[3]], net_size)
                zzz = np.clip(zoom(zzz, dilations), 0, 1)
                # If we have some rounding problems from zoom, pad with 0 for the rounding errors.
                if zzz.shape < croppedIntensityMask.shape:
                    zzz = scale(zzz, np.array([zzz.shape[0], zzz.shape[1], zzz.shape[2]]),
                                  np.array([croppedIntensityMask.shape[0], croppedIntensityMask.shape[1], croppedIntensityMask.shape[2]]))
                if zzz.shape > croppedIntensityMask.shape:
                    zzz = crop(zzz, np.array([croppedIntensityMask.shape[0], croppedIntensityMask.shape[1], croppedIntensityMask.shape[2]]),
                                np.array([zzz.shape[0], zzz.shape[1], zzz.shape[2]]))
            else:
                zzz = np.squeeze(np.squeeze(sess.run(net, {lesionPlaceholder: expandedLesion}), 4), 0)
                zzz = np.clip(zoom(zzz, dilations), 0, 1)
            lesionPriorVAE = np.pad((zzz * croppedIntensityMask), paddings, mode='constant')[ mask ]

            if hasattr( visualizer_, 'show_flag' ): 
                tmp = np.zeros( imageSize ); tmp[ mask ] = lesionPriorVAE
                visualizer_.show( probabilities=tmp, title="Lesion prior using VAE only", window_id="Lesion prior using VAE only" )


            # Sample from the mean and variance, conditioned on the data and the lesion segmentation 
            bestMean, bestVariance, _ = fitGMMParameters( data, 
                                                          lesion[ mask ].reshape( -1, 1 ), 
                                                          [ 1 ], useDiagonalCovarianceMatrices,
                                                          hyperMeans=np.array( [ hyperMean ] ), 
                                                          hyperMeansNumberOfMeasurements=np.array( [ hyperMeanNumberOfMeasurements ] ),
                                                          hyperVariances=np.array( [ hyperVariance ] ), 
                                                          hyperVariancesNumberOfMeasurements=np.array( [ hyperVarianceNumberOfMeasurements ] )
                                                        )
            bestMean, bestVariance = bestMean[ 0 ], bestVariance[ 0 ]
            N = lesion[ mask ].sum()
            
            # Murphy, page 134 with v0 = hyperVarianceNumberOfMeasurements - numberOfContrasts - 1
            variance = invwishart.rvs( N + hyperVarianceNumberOfMeasurements - numberOfContrasts - 1, 
                                      bestVariance * ( hyperVarianceNumberOfMeasurements + N + 1 ) )
            if useDiagonalCovarianceMatrices:
                variance = np.diag( np.diag( variance ) )
            mean = np.random.multivariate_normal( bestMean, variance / ( hyperMeanNumberOfMeasurements + N ) ).reshape( -1, 1 )
            
            
            
            
            # Sample from the lesion segmentation, conditioned on the data and the VAE latent variables
            # (Implementation-wise the latter is encoded in the VAE prior). At the same time we also
            # compute the full posterior of each structure, which is at the end the thing we're averaging
            # over (i.e., the reason why we're sampling)
            effectivePriors = np.array( priors / 65535, dtype=np.float32 )
            effectivePriors[:, lesionStructureNumber ] *= lesionPriorVAE
            if hasattr( visualizer_, 'show_flag' ): 
                tmp = np.zeros( imageSize ); tmp[ mask ] = effectivePriors[:, lesionStructureNumber ]
                visualizer_.show( probabilities=tmp, title="Lesion prior using VAE and atlas together",
                                window_id="Lesion prior using VAE and atlas together" )
            if False:
                # This is what Stefano had in his code. Hard to see what sort of model this might correspond to
                normalizer = np.sum(effectivePriors, axis=1) + eps
                effectivePriors = effectivePriors / np.expand_dims(normalizer, 1)
            else:
                # Generative model where the atlas generates *candidate* lesions, and the VAE prior is sampled 
                # from *only within the candidates*.
                numberOfStructures = priors.shape[-1]
                otherStructureNumbers = [ i for i in range( numberOfStructures ) if i != lesionStructureNumber ]
                otherStructurePriors = effectivePriors[ :, otherStructureNumbers ]
                otherStructurePriors /= ( np.sum( otherStructurePriors, axis=1 ).reshape( -1, 1 ) + eps )
                effectivePriors[ :, otherStructureNumbers ] = otherStructurePriors * \
                                                            ( 1 - effectivePriors[:, lesionStructureNumber ].reshape( -1, 1 ) )
            likelihoods[ :, lesionStructureNumber ] = getGaussianLikelihoods( data, mean, variance )
            posteriors = effectivePriors * likelihoods
            posteriors /= np.expand_dims( np.sum( posteriors, axis=1 ) + eps, 1 )
            sample = np.random.rand( numberOfVoxels ) <= posteriors[ :, lesionStructureNumber ] 
            lesion = np.zeros( imageSize ); lesion[ mask ] = sample

            visualizer_.show( image_list=[ lesion ], title="Lesion sample", window_id="Lesion sample" )
            
            
        # Collect data after burn in steps
        if sweepNumber >= numberOfBurnInSteps_:
            print('Sample ' + str( sweepNumber + 1 - numberOfBurnInSteps_ ) + ' times')
            averagePosteriors += posteriors / numberOfSamplingSteps_
        else:
            print('Burn-in ' + str( sweepNumber + 1) + ' times')


    #
    visualizer_.show_movie( window_id="Lesion prior using VAE only" )
    visualizer_.show_movie( window_id="Lesion sample" )


    # Return
    return averagePosteriors



def samsegmentLesion( imageFileNames, atlasDir, savePath,
                      transformedTemplateFileName=None, 
                      userModelSpecifications={}, userOptimizationOptions={},
                      lesionSearchString='Lesion',
                      numberOfSamplingSteps=50, numberOfBurnInSteps=50,
                      visualizer=None, saveHistory=False, saveMesh=False,
                      targetIntensity=None, targetSearchStrings=None,
                      intensityMaskingPattern=None, intensityMaskingSearchString=None,
                      lesionThreshold=None ): 
  
    # We can't pass on user-specified options to plugins, so let's do it using global variables
    # (using underscore notation to stress they're global variables )
    global hyperMeans_, hyperMeansNumberOfMeasurements_, hyperVariances_, hyperVariancesNumberOfMeasurements_, \
           atlasDir_, lesionClassNumber_ , intensityMaskingClassNumber_, numberOfSamplingSteps_, \
           numberOfBurnInSteps_, intensityMaskingPattern_, visualizer_
    
    atlasDir_ = atlasDir
    lesionClassNumber_ = getClassNumber( lesionSearchString )
    intensityMaskingClassNumber_ = getClassNumber( intensityMaskingSearchString )
    numberOfSamplingSteps_ = numberOfSamplingSteps
    numberOfBurnInSteps_ = numberOfBurnInSteps
    intensityMaskingPattern_ = intensityMaskingPattern
    visualizer_ = visualizer


    # Now call samsegment with plugins
    samsegment( imageFileNames, atlasDir, savePath,
                userModelSpecifications=userModelSpecifications, 
                userOptimizationOptions=userOptimizationOptions,
                visualizer=visualizer,
                targetIntensity=targetIntensity, 
                targetSearchStrings=targetSearchStrings,
                hyperpriorPlugin=defineHyperparameters,
                posteriorPlugin=getMCMCPosteriors,
                posteriorPluginVariables=[ 'modelSpecifications.useDiagonalCovarianceMatrices', 
                                            'mask', 'voxelSpacing', 'transform' ],
                threshold=lesionThreshold, thresholdSearchString=lesionSearchString )

