import os
import math
import numpy as np
import freesurfer as fs
import freesurfer.gems as gems
from freesurfer.gems import kvlReadCompressionLookupTable, kvlReadSharedGMMParameters
from freesurfer.samseg.figures import initVisualizer
from freesurfer.samseg.utilities import requireNumpyArray, Specification
from freesurfer.samseg.bias_correction import projectKroneckerProductBasisFunctions, backprojectKroneckerProductBasisFunctions, \
    computePrecisionOfKroneckerProductBasisFunctions
from freesurfer.samseg.register_atlas import registerAtlas
import logging
import pickle


#eps = np.spacing(1)
eps = np.finfo( float ).eps


def showImage( data ):
    range = ( data.min(), data.max() )

    Nx = data.shape[0]
    Ny = data.shape[1]
    Nz = data.shape[2]

    x = round( Nx / 2 )
    y = round( Ny / 2 )
    z = round( Nz / 2 )

    xySlice = data[ :, :, z ]
    xzSlice = data[ :, y, : ]
    yzSlice = data[ x, :, : ]

    patchedSlices = np.block( [ [ xySlice, xzSlice], [ yzSlice.T, np.zeros( ( Nz, Nz ) ) + range[ 0 ] ] ] )

    import matplotlib.pyplot as plt   # avoid importing matplotlib by default
    plt.imshow( patchedSlices.T, cmap=plt.cm.gray, vmin=range[ 0 ], vmax=range[ 1 ] )
    #plt.gray()
    #plt.imshow( patchedSlices.T, vmin=range[ 0 ], vmax=range[ 1 ] )
    #plt.show()
    plt.axis( 'off' )



def getModelSpecifications( atlasDir, userModelSpecifications={} ):
  
    # Create default model specifications as a dictionary
    FreeSurferLabels, names, colors = kvlReadCompressionLookupTable( os.path.join( atlasDir, 'compressionLookupTable.txt') )
    sharedGMMParameters = kvlReadSharedGMMParameters( os.path.join( atlasDir, 'sharedGMMParameters.txt') )

    modelSpecifications = {
        'FreeSurferLabels': FreeSurferLabels,
        'atlasFileName': os.path.join(atlasDir, 'atlas_level2.txt.gz'),
        'names': names,
        'colors': colors,
        'sharedGMMParameters': sharedGMMParameters,
        'useDiagonalCovarianceMatrices': True,
        'brainMaskingSmoothingSigma': 3.0,  # sqrt of the variance of a Gaussian blurring kernel
        'brainMaskingThreshold': 0.01,
        'K': 0.1,  # stiffness of the mesh
        'biasFieldSmoothingKernelSize': 50,  # distance in mm of sinc function center to first zero crossing
    }

    modelSpecifications.update( userModelSpecifications )

    return modelSpecifications


def getOptimizationOptions( atlasDir, userOptimizationOptions={} ):
  
    # Create default optimization options as a dictionary
    optimizationOptions = {
        'maximumNumberOfDeformationIterations': 20,
        'absoluteCostPerVoxelDecreaseStopCriterion': 1e-4,
        'verbose': False,
        'maximalDeformationStopCriterion': 0.001,  # measured in pixels
        'lineSearchMaximalDeformationIntervalStopCriterion': 0.001,
        'maximalDeformationAppliedStopCriterion': 0.0,
        'BFGSMaximumMemoryLength': 12,
        'multiResolutionSpecification': 
          [ 
            { 'atlasFileName': os.path.join( atlasDir, 'atlas_level1.txt.gz' ),
              'targetDownsampledVoxelSpacing': 2.0,
              'maximumNumberOfIterations': 100,
              'estimateBiasField': True 
            },
            { 'atlasFileName': os.path.join( atlasDir, 'atlas_level2.txt.gz' ),
              'targetDownsampledVoxelSpacing': 1.0,
              'maximumNumberOfIterations': 100,
              'estimateBiasField': True 
            }
          ]
    }
            
            
    # Over-write with any uper specified options. The 'multiResolutionSpecification' key has as value a list
    # of dictionaries which we shouldn't just over-write, but rather update themselves, so this is special case
    userOptimizationOptionsCopy = userOptimizationOptions.copy()
    key = 'multiResolutionSpecification'
    if key in userOptimizationOptionsCopy:
        userList = userOptimizationOptionsCopy[ key ]
        defaultList = optimizationOptions[ key ]
        for levelNumber in range( len( defaultList ) ):
            if levelNumber < len( userList ):
                defaultList[ levelNumber ].update( userList[ levelNumber ] )
            else:
                del defaultList[ levelNumber ]
        del userOptimizationOptionsCopy[ key ]
    optimizationOptions.update( userOptimizationOptionsCopy )  
        
    return optimizationOptions

 
  
def readCroppedImages( imageFileNames, transformedTemplateFileName ):

  # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
  # translation, rotation, scaling, and skewing) as well - this transformation will later be used
  # to initially transform the location of the atlas mesh's nodes into the coordinate system of the image.
  imageBuffers = []
  for imageFileName in imageFileNames:
      # Get the pointers to image and the corresponding transform
      image = gems.KvlImage( imageFileName, transformedTemplateFileName )
      transform = image.transform_matrix
      cropping = image.crop_slices
      imageBuffers.append( image.getImageBuffer() )

  imageBuffers = np.transpose( imageBuffers, axes=[1, 2, 3, 0] )


  # Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels, 
  # downsampling steps etc in mm.
  nonCroppedImage = gems.KvlImage( imageFileNames[0] )
  imageToWorldTransformMatrix = nonCroppedImage.transform_matrix.as_numpy_array
  voxelSpacing = np.sum( imageToWorldTransformMatrix[ 0:3, 0:3 ] ** 2, axis=0 ) ** ( 1/2 )

  #
  return imageBuffers, transform, voxelSpacing, cropping


def maskOutBackground( imageBuffers, atlasFileName, transform,
                       brainMaskingSmoothingSigma, brainMaskingThreshold, 
                       visualizer=None, maskOutZeroIntensities=True ):
    # Setup a null visualizer if necessary
    if visualizer is None: visualizer = initVisualizer(False, False)
  
    # Read the affinely coregistered atlas mesh (in reference position)
    mesh = getMesh( atlasFileName, transform )
  
    # Mask away uninteresting voxels. This is done by a poor man's implementation of a dilation operation on
    # a non-background class mask; followed by a cropping to the area covered by the mesh (needed because
    # otherwise there will be voxels in the data with prior probability zero of belonging to any class)
    imageSize = imageBuffers.shape[ 0:3 ]
    labelNumber = 0
    backgroundPrior = mesh.rasterize_1a( imageSize, labelNumber )

    # Threshold background prior at 0.5 - this helps for atlases built from imperfect (i.e., automatic)
    # segmentations, whereas background areas don't have zero probability for non-background structures
    backGroundThreshold = 2 ** 8
    backGroundPeak = 2 ** 16 - 1
    backgroundPrior = np.ma.filled( np.ma.masked_greater( backgroundPrior, backGroundThreshold ),
                                    backGroundPeak ).astype( np.float32 )

    visualizer.show( probabilities=backgroundPrior, images=imageBuffers, window_id='samsegment background', 
                     title='Background Priors' )
    
    smoothingSigmas = [1.0 * brainMaskingSmoothingSigma] * 3
    smoothedBackgroundPrior = gems.KvlImage.smooth_image_buffer( backgroundPrior, smoothingSigmas )
    visualizer.show( probabilities=smoothedBackgroundPrior, window_id='samsegment smoothed', 
                     title='Smoothed Background Priors')


    # 65535 = 2^16 - 1. priors are stored as 16bit ints
    # To put the threshold in perspective: for Gaussian smoothing with a 3D isotropic kernel with variance
    # diag( sigma^2, sigma^2, sigma^2 ) a single binary "on" voxel at distance sigma results in a value of
    # 1/( sqrt(2*pi)*sigma )^3 * exp( -1/2 ).
    # More generally, a single binary "on" voxel at some Eucledian distance d results in a value of
    # 1/( sqrt(2*pi)*sigma )^3 * exp( -1/2*d^2/sigma^2 ). Turning this around, if we threshold this at some
    # value "t", a single binary "on" voxel will cause every voxel within Eucledian distance
    #
    #   d = sqrt( -2*log( t * ( sqrt(2*pi)*sigma )^3 ) * sigma^2 )
    #
    # of it to be included in the mask.
    #
    # As an example, for 1mm isotropic data, the choice of sigma=3 and t=0.01 yields ... complex value ->
    # actually a single "on" voxel will then not make any voxel survive, as the normalizing constant (achieved
    # at Mahalanobis distance zero) is already < 0.01
    brainMaskThreshold = 65535.0 * ( 1.0 - brainMaskingThreshold )
    brainMask = np.ma.less( smoothedBackgroundPrior, brainMaskThreshold )

    # Crop to area covered by the mesh
    alphas = mesh.alphas
    areaCoveredAlphas = [[0.0, 1.0]] * alphas.shape[0]
    mesh.alphas = areaCoveredAlphas  # temporary replacement of alphas
    areaCoveredByMesh = mesh.rasterize_1b( imageSize, 1 )
    mesh.alphas = alphas  # restore alphas
    brainMask = np.logical_and( brainMask, areaCoveredByMesh )

    # If a pixel has a zero intensity in any of the contrasts, that is also masked out across all contrasts
    if maskOutZeroIntensities:
        numberOfContrasts = imageBuffers.shape[-1]
        for contrastNumber in range( numberOfContrasts ):
            brainMask *= imageBuffers[ :, :, :, contrastNumber ] > 0
            
    # Mask the images
    maskedImageBuffers = imageBuffers.copy()
    maskedImageBuffers[ np.logical_not( brainMask ), : ] = 0

    #
    return maskedImageBuffers, brainMask


def getBiasFieldBasisFunctions( imageSize, smoothingKernelSize ):
  
    # Our bias model is a linear combination of a set of basis functions. We are using so-called
    # "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
    # Transform.
    biasFieldBasisFunctions = []
    for dimensionNumber in range( 3 ):
        N = imageSize[ dimensionNumber ]
        delta = smoothingKernelSize[ dimensionNumber ]
        M = math.ceil( N / delta ) + 1
        Nvirtual = ( M - 1 ) * delta
        js = [ (index + 0.5) * math.pi / Nvirtual for index in range( N ) ]
        scaling = [ math.sqrt( 2 / Nvirtual ) ] * M
        scaling[ 0 ] /= math.sqrt( 2 )
        A = np.array( [ [ math.cos( freq * m ) * scaling[ m ] for m in range( M ) ] for freq in js ] )
        biasFieldBasisFunctions.append( A )
        
            
    return biasFieldBasisFunctions



def getMesh( meshCollectionFileName, 
             transform=None, 
             K=None, 
             initialDeformation=None, initialDeformationMeshCollectionFileName=None, 
             returnInitialDeformationApplied=False ):

    # Get the mesh
    mesh_collection = gems.KvlMeshCollection()
    mesh_collection.read( meshCollectionFileName )
    if K:
        mesh_collection.k = K
    if transform:
        mesh_collection.transform( transform )
    else:
        transform = gems.KvlTransform( requireNumpyArray( np.eye( 4 ) ) )
    mesh = mesh_collection.reference_mesh
    
    
    # See if we need to warp it
    estimatedNodeDeformationInTemplateSpace = None
    if initialDeformation is not None:
        #
        if initialDeformationMeshCollectionFileName is None:
            initialDeformationMeshCollectionFileName = meshCollectionFileName
      
        # Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
        nodePositions = mesh.points
        nodePositionsInTemplateSpace = mapPositionsFromSubjectToTemplateSpace( nodePositions, transform )

        # Get the estimated warp in template space
        [ estimatedNodeDeformationInTemplateSpace, estimated_averageDistance, estimated_maximumDistance ] = gems.kvlWarpMesh(
            initialDeformationMeshCollectionFileName,
            initialDeformation,
            meshCollectionFileName
        )

        # Apply this warp on the mesh node positions in template space, and transform into current space
        desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace
        desiredNodePositions = mapPositionsFromTemplateToSubjectSpace( desiredNodePositionsInTemplateSpace, transform )
        mesh.points = requireNumpyArray( desiredNodePositions )

    # Return what we got
    if returnInitialDeformationApplied:
        if estimatedNodeDeformationInTemplateSpace is None:
            estimatedNodeDeformationInTemplateSpace = np.zeros_like( mesh.points )
        return mesh, estimatedNodeDeformationInTemplateSpace
    else:
        return mesh


def getDownSampledModel( imageBuffers, mask, atlasFileName, K, transform, biasFieldBasisFunctions, 
                         downSamplingFactors, initialDeformation=None, initialDeformationMeshCollectionFileName=None ):
  
    # Downsample the images and basis functions
    numberOfContrasts = imageBuffers.shape[-1]
    downSampledMask = mask[ ::downSamplingFactors[0], ::downSamplingFactors[1], ::downSamplingFactors[2] ]
    downSampledImageBuffers = np.zeros( downSampledMask.shape + (numberOfContrasts,), order='F' )
    for contrastNumber in range( numberOfContrasts ):
        #logger.debug('first time contrastNumber=%d', contrastNumber)
        downSampledImageBuffers[ :, :, :, contrastNumber ] = imageBuffers[ ::downSamplingFactors[0],
                                                                           ::downSamplingFactors[1],
                                                                           ::downSamplingFactors[2],
                                                                           contrastNumber ]

    downSampledBiasFieldBasisFunctions = [ np.array( biasFieldBasisFunction[ ::downSamplingFactor ] )
                                                  for biasFieldBasisFunction, downSamplingFactor in
                                                  zip( biasFieldBasisFunctions, downSamplingFactors ) ]
    
    # Compute the resulting transform, taking into account the downsampling
    downSamplingTransformMatrix = np.diag( 1. / downSamplingFactors )
    downSamplingTransformMatrix = np.pad( downSamplingTransformMatrix, (0, 1), mode='constant', constant_values=0 )
    downSamplingTransformMatrix[3][3] = 1
    downSampledTransform = gems.KvlTransform( requireNumpyArray( downSamplingTransformMatrix @ transform.as_numpy_array ) )

    # Get the mesh
    downSampledMesh, downSampledInitialDeformationApplied = getMesh( atlasFileName, downSampledTransform, K, 
                                                                     initialDeformation, initialDeformationMeshCollectionFileName, 
                                                                     returnInitialDeformationApplied=True ) 
  
    return downSampledImageBuffers, downSampledMask, downSampledMesh, downSampledInitialDeformationApplied, \
           downSampledTransform, downSampledBiasFieldBasisFunctions
         
         
def mapPositionsFromSubjectToTemplateSpace( positions, transform ):

    #
    tmp = np.linalg.solve( transform.as_numpy_array,
                           np.pad( positions, ( (0, 0), (0, 1) ), mode='constant', constant_values=1 ).T ).T
    return tmp[:, 0:3]

def mapPositionsFromTemplateToSubjectSpace( positions, transform ):
  
    #
    tmp = ( transform.as_numpy_array @ \
            np.pad( positions, ((0, 0), (0, 1)), 'constant', constant_values=1 ).T ).T
    return tmp[:, 0:3]


def getBiasFields( biasFieldCoefficients, biasFieldBasisFunctions,  mask=None ):

    #
    numberOfContrasts = biasFieldCoefficients.shape[-1]
    imageSize = tuple( [ functions.shape[0] for functions in biasFieldBasisFunctions ] )
    biasFields = np.zeros( imageSize + (numberOfContrasts,), order='F' )
    for contrastNumber in range( numberOfContrasts ):
        biasField = backprojectKroneckerProductBasisFunctions(
              biasFieldBasisFunctions, biasFieldCoefficients[ :, contrastNumber ] )
        if mask is not None:
            biasField *= mask
        biasFields[ :, :, :, contrastNumber ] = biasField

    return biasFields


def getGaussianLikelihoods( data, mean, variance ):
  
  #
  numberOfContrasts = data.shape[1]
  
  L = np.linalg.cholesky( variance )
  tmp = np.linalg.solve( L, data.T - mean )
  squaredMahalanobisDistances = np.sum( tmp ** 2, axis=0 )
  sqrtDeterminantOfVariance = np.prod( np.diag( L ) )
  scaling = 1.0 / (2 * np.pi) ** ( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance
  gaussianLikelihoods = np.exp( squaredMahalanobisDistances * -0.5 ) * scaling
  return gaussianLikelihoods.T


def getGaussianPosteriors( data, classPriors, means, variances, mixtureWeights, numberOfGaussiansPerClass ):

    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    numberOfContrasts = data.shape[1]
    numberOfVoxels = data.shape[0]
    
    gaussianPosteriors = np.zeros( ( numberOfVoxels, numberOfGaussians ), order='F' )
    for classNumber in range( numberOfClasses ):
        classPrior = classPriors[ :, classNumber ]
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
        for componentNumber in range( numberOfComponents ):
            gaussianNumber = sum( numberOfGaussiansPerClass[ :classNumber ] ) + componentNumber
            mean = np.expand_dims( means[ gaussianNumber, : ], 1 )
            variance = variances[ gaussianNumber, :, : ]

            gaussianLikelihoods = getGaussianLikelihoods( data, mean, variance )
            gaussianPosteriors[ :, gaussianNumber ] = gaussianLikelihoods * ( mixtureWeights[ gaussianNumber ] * classPrior )
    normalizer = np.sum( gaussianPosteriors, axis=1 ) + eps
    gaussianPosteriors = gaussianPosteriors / np.expand_dims( normalizer, 1 )

    minLogLikelihood = -np.sum( np.log( normalizer ) )
 
    return gaussianPosteriors, minLogLikelihood
  
  
def getLikelihoods( data, means, variances, mixtureWeights, numberOfGaussiansPerClass, fractionsTable ):
    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    numberOfContrasts = data.shape[1]
    numberOfVoxels = data.shape[0]
    numberOfStructures = fractionsTable.shape[1]

    #   
    likelihoods = np.zeros( ( numberOfVoxels, numberOfStructures ), dtype=np.float64 )
    for classNumber in range( numberOfClasses ):
      
        # Compute likelihood for this class
        classLikelihoods = np.zeros( numberOfVoxels )
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
        for componentNumber in range( numberOfComponents ):
            gaussianNumber = sum( numberOfGaussiansPerClass[ :classNumber ] ) + componentNumber
            mean = np.expand_dims( means[ gaussianNumber, : ], 1 )
            variance = variances[ gaussianNumber, :, : ]
            mixtureWeight = mixtureWeights[ gaussianNumber ]
            
            gaussianLikelihoods = getGaussianLikelihoods( data, mean, variance )
            classLikelihoods += gaussianLikelihoods * mixtureWeight
            
        # Add contribution to the actual structures            
        for structureNumber in range( numberOfStructures ):
            fraction = fractionsTable[ classNumber, structureNumber ]
            if fraction < 1e-10: 
                continue
            likelihoods[ :, structureNumber ] += classLikelihoods * fraction
  
    #
    return likelihoods
  
  
def getPosteriors( data, priors, means, variances, mixtureWeights, numberOfGaussiansPerClass, fractionsTable ):

    # Weight likelihood against prior and normalize
    posteriors = getLikelihoods( data, means, variances, mixtureWeights, 
                                 numberOfGaussiansPerClass, fractionsTable ) * priors
    normalizer = np.sum( posteriors, axis=1 ) + eps
    posteriors = posteriors / np.expand_dims( normalizer, 1 )

    return posteriors
  
  
  
  

def getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts,
                            hyperMeans=None, hyperMeansNumberOfMeasurements=None, 
                            hyperVariances=None, hyperVariancesNumberOfMeasurements=None,
                            hyperMixtureWeights=None, hyperMixtureWeightsNumberOfMeasurements=None,
                          ):
  
    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )

  
    # Make sure the hyperparameters are defined
    if hyperMeans is None:
        hyperMeans = np.zeros( ( numberOfGaussians, numberOfContrasts ) )
    if hyperMeansNumberOfMeasurements is None:
        hyperMeansNumberOfMeasurements = np.zeros( numberOfGaussians )
    else:
        hyperMeansNumberOfMeasurements = hyperMeansNumberOfMeasurements.copy()
    if hyperVariances is None:
        hyperVariances = np.tile( np.eye( numberOfContrasts ), ( numberOfGaussians, 1, 1 ) )
    if hyperVariancesNumberOfMeasurements is None:
        hyperVariancesNumberOfMeasurements = np.zeros( numberOfGaussians )
    else:
        hyperVariancesNumberOfMeasurements = hyperVariancesNumberOfMeasurements.copy()
    if hyperMixtureWeights is None:
        hyperMixtureWeights = np.ones( numberOfGaussians )
        for classNumber in range( numberOfClasses ):
            # mixture weights are normalized (those belonging to one mixture sum to one)
            numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
            gaussianNumbers = np.array( np.sum( numberOfGaussiansPerClass[ :classNumber ] ) + \
                                        np.array( range( numberOfComponents ) ), dtype=np.uint32 )
            hyperMixtureWeights[ gaussianNumbers ] /= np.sum( hyperMixtureWeights[ gaussianNumbers ] )
    if hyperMixtureWeightsNumberOfMeasurements is None:
        hyperMixtureWeightsNumberOfMeasurements = np.zeros( numberOfClasses )
    else:
        hyperMixtureWeightsNumberOfMeasurements = hyperMixtureWeightsNumberOfMeasurements.copy()

    # Making sure the inverse-Wishart is normalizable (flat or peaked around hyperVariances)
    # requires that hyperVarianceNumberOfMeasurements is not smaller than (numberOfContrasts-1)
    # for any Gaussian. However, in order to prevent numerical errors with near-zero variances
    # (which can happen when almost no voxels are associated with a Gaussian in the EM algorithm,
    # due to e.g., tiny mixture weight), we use (numberOfContrasts-1)+1 instead.
    threshold = ( numberOfContrasts - 1 ) + 1 + eps
    hyperVariancesNumberOfMeasurements[ hyperVariancesNumberOfMeasurements < threshold ] = threshold

    if False:
        print( 'hyperMeans: ', hyperMeans )
        print( 'hyperMeansNumberOfMeasurements: ', hyperMeansNumberOfMeasurements )
        print( 'hyperVariances: ', hyperVariances )
        print( 'hyperVariancesNumberOfMeasurements: ', hyperVariancesNumberOfMeasurements )
        print( 'hyperMixtureWeights: ', hyperMixtureWeights )
        print( 'hyperMixtureWeightsNumberOfMeasurements: ', hyperMixtureWeightsNumberOfMeasurements )


    #
    return hyperMeans, hyperMeansNumberOfMeasurements, \
           hyperVariances, hyperVariancesNumberOfMeasurements, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements


def getDownsampledHyperparameters( hyperMeans, hyperMeansNumberOfMeasurements, 
                                   hyperVariances, hyperVariancesNumberOfMeasurements,
                                   hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements,
                                   downSamplingFactors ):
  
    return hyperMeans, \
           hyperMeansNumberOfMeasurements / np.prod( downSamplingFactors ) \
              if hyperMeansNumberOfMeasurements is not None else None, \
           hyperVariances, \
           hyperVariancesNumberOfMeasurements / np.prod( downSamplingFactors ) \
              if hyperVariancesNumberOfMeasurements is not None else None, \
           hyperMixtureWeights, \
           hyperMixtureWeightsNumberOfMeasurements  / np.prod( downSamplingFactors ) \
              if hyperMixtureWeightsNumberOfMeasurements is not None else None



def fitGMMParameters( data, gaussianPosteriors, numberOfGaussiansPerClass, useDiagonalCovarianceMatrices,
                      hyperMeans=None, hyperMeansNumberOfMeasurements=None, 
                      hyperVariances=None, hyperVariancesNumberOfMeasurements=None,
                      hyperMixtureWeights=None, hyperMixtureWeightsNumberOfMeasurements=None 
                    ):
  
    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    numberOfContrasts = data.shape[1]
    
    # Make sure the hyperparameters are defined and valid
    hyperMeans, hyperMeansNumberOfMeasurements, \
           hyperVariances, hyperVariancesNumberOfMeasurements, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements = \
               getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts,
                                        hyperMeans, hyperMeansNumberOfMeasurements, 
                                        hyperVariances, hyperVariancesNumberOfMeasurements,
                                        hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements )

    # Means and variances
    means = np.zeros( ( numberOfGaussians, numberOfContrasts ) )
    variances = np.zeros( ( numberOfGaussians, numberOfContrasts, numberOfContrasts ) )
    for gaussianNumber in range( numberOfGaussians ):
        posterior = gaussianPosteriors[ :, gaussianNumber ].reshape( -1, 1 )
        hyperMean = np.expand_dims( hyperMeans[ gaussianNumber, : ], 1 )
        hyperMeanNumberOfMeasurements = hyperMeansNumberOfMeasurements[ gaussianNumber ]
        hyperVariance = hyperVariances[ gaussianNumber, :, : ]
        hyperVarianceNumberOfMeasurements = hyperVariancesNumberOfMeasurements[ gaussianNumber ]
        
        mean = ( data.T @ posterior + hyperMean * hyperMeanNumberOfMeasurements ) \
               / ( np.sum( posterior ) + hyperMeanNumberOfMeasurements )
        tmp = data - mean.T
        variance = ( tmp.T @ (tmp * posterior) + \
                    hyperMeanNumberOfMeasurements * ( ( mean - hyperMean ) @ ( mean - hyperMean ).T ) + \
                    hyperVariance * hyperVarianceNumberOfMeasurements ) \
                    / ( np.sum(posterior) + 1 + hyperVarianceNumberOfMeasurements )
        if useDiagonalCovarianceMatrices:
            # Force diagonal covariance matrices
            variance = np.diag( np.diag( variance ) )
        variances[ gaussianNumber, :, : ] = variance
        means[ gaussianNumber, : ] = mean.T
        
    # Mixture weights    
    mixtureWeights = np.sum( gaussianPosteriors + eps, axis=0 )
    for classNumber in range( numberOfClasses ):
        # mixture weights are normalized (those belonging to one mixture sum to one)
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
        gaussianNumbers = np.array( np.sum( numberOfGaussiansPerClass[ :classNumber ] ) + \
                                    np.array( range( numberOfComponents ) ), dtype=np.uint32 )
        
        mixtureWeights[ gaussianNumbers ] += hyperMixtureWeights[ gaussianNumbers ] * \
                                             hyperMixtureWeightsNumberOfMeasurements[ classNumber ]
        mixtureWeights[ gaussianNumbers ] /= np.sum( mixtureWeights[ gaussianNumbers ] )

    #
    return means, variances, mixtureWeights


  
  
def evaluateMinLogPriorOfGMMParameters( means, variances, mixtureWeights, 
                                        numberOfGaussiansPerClass, 
                                        hyperMeans=None, hyperMeansNumberOfMeasurements=None, 
                                        hyperVariances=None, hyperVariancesNumberOfMeasurements=None,
                                        hyperMixtureWeights=None, hyperMixtureWeightsNumberOfMeasurements=None 
                                       ):
  
    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    numberOfContrasts = means.shape[1]
    
    # Make sure the hyperparameters are defined and valid
    hyperMeans, hyperMeansNumberOfMeasurements, \
           hyperVariances, hyperVariancesNumberOfMeasurements, \
           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements = \
               getFullHyperparameters( numberOfGaussiansPerClass, numberOfContrasts,
                                        hyperMeans, hyperMeansNumberOfMeasurements, 
                                        hyperVariances, hyperVariancesNumberOfMeasurements,
                                        hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements )

    # 
    minLogPrior = 0
    for gaussianNumber in range( numberOfGaussians ):
        mean = np.expand_dims( means[ gaussianNumber, : ], 1 )
        variance = variances[ gaussianNumber, :, : ]

        hyperMean = np.expand_dims( hyperMeans[ gaussianNumber, : ], 1 )
        hyperMeanNumberOfMeasurements = hyperMeansNumberOfMeasurements[ gaussianNumber ]
        hyperVariance = hyperVariances[ gaussianNumber, :, : ]
        hyperVarianceNumberOfMeasurements = hyperVariancesNumberOfMeasurements[ gaussianNumber ]

        # -log N( mean | hyperMean, variance / hyperMeanNumberOfMeasurements )
        L = np.linalg.cholesky( variance ) # variance = L @ L.T
        halfOfLogDetVariance = np.sum( np.log( np.diag( L ) ) )
        tmp = np.linalg.solve( L, mean - hyperMean )
        squaredMahalanobisDistance = np.sum( tmp * tmp )
        minLogPrior += squaredMahalanobisDistance * hyperMeanNumberOfMeasurements / 2 + halfOfLogDetVariance

        # -log IW( variance | hyperVariance * hyperVarianceNumberOfMeasurements, 
        #                     hyperVarianceNumberOfMeasurements - numberOfContrasts - 1 )
        #
        hyperL = np.linalg.cholesky( hyperVariance ) # hyperVariance = hyperL @ hyperL.T
        halfOfLogDetHyperVariance = np.sum( np.log( np.diag( hyperL ) ) )
        tmp = np.linalg.solve( L, hyperL )
        minLogPrior += np.trace( tmp @ tmp.T ) * hyperVarianceNumberOfMeasurements / 2 + \
                       hyperVarianceNumberOfMeasurements * halfOfLogDetVariance - \
                       ( hyperVarianceNumberOfMeasurements - numberOfContrasts - 1 ) * halfOfLogDetHyperVariance


    for classNumber in range( numberOfClasses ):
        # -log Dir( weights | hyperMixtureWeights * hyperMixtureWeightNumberOfMeasurements + 1 )
        hyperMixtureWeightNumberOfMeasurements = hyperMixtureWeightsNumberOfMeasurements[ classNumber ]
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
        for componentNumber in range( numberOfComponents ):
            gaussianNumber = sum( numberOfGaussiansPerClass[ :classNumber ] ) + componentNumber
            mixtureWeight = mixtureWeights[ gaussianNumber ]
            hyperMixtureWeight = hyperMixtureWeights[ gaussianNumber ]

            # minLogPrior -= hyperMixtureWeight * hyperMixtureWeightNumberOfMeasurements * np.log( mixtureWeight )
            #
            # I'm using Stirling's approximation on the normalizing constant (beta function) just in the same way
            # as in Appendix C of Van Leemput TMI 2009
            minLogPrior += hyperMixtureWeightNumberOfMeasurements * \
                           hyperMixtureWeight * ( np.log( hyperMixtureWeight + eps ) - np.log( mixtureWeight + eps ) )
            
    #
    return minLogPrior
  
  
  

def fitBiasFieldParameters( imageBuffers, gaussianPosteriors, means, variances, biasFieldBasisFunctions, mask ):

    # Bias field correction: implements Eq. 8 in the paper
    #    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999

    #
    numberOfGaussians = means.shape[ 0 ]
    numberOfContrasts = means.shape[ 1 ]
    numberOfBasisFunctions = [ functions.shape[1] for functions in biasFieldBasisFunctions ]
    numberOf3DBasisFunctions = np.prod( numberOfBasisFunctions )

    # Set up the linear system lhs * x = rhs
    precisions = np.zeros_like( variances )
    for gaussianNumber in range( numberOfGaussians ):
        precisions[ gaussianNumber, :, : ] = np.linalg.inv( variances[ gaussianNumber, :, : ] ).reshape(
            ( 1, numberOfContrasts, numberOfContrasts ) )
        
    lhs = np.zeros( ( numberOf3DBasisFunctions * numberOfContrasts, 
                      numberOf3DBasisFunctions * numberOfContrasts ) )  # left-hand side of linear system
    rhs = np.zeros( ( numberOf3DBasisFunctions * numberOfContrasts, 1 ) )  # right-hand side of linear system
    weightsImageBuffer = np.zeros( mask.shape )
    tmpImageBuffer = np.zeros( mask.shape )
    for contrastNumber1 in range( numberOfContrasts ):
        #logger.debug('third time contrastNumber=%d', contrastNumber)
        contrast1Indices = np.arange( 0, numberOf3DBasisFunctions ) + \
                              contrastNumber1 * numberOf3DBasisFunctions
        
        tmp = np.zeros( gaussianPosteriors.shape[0] )
        for contrastNumber2 in range( numberOfContrasts ):
            contrast2Indices = np.arange( 0, numberOf3DBasisFunctions ) + \
                                contrastNumber2 * numberOf3DBasisFunctions
          
            classSpecificWeights = gaussianPosteriors * precisions[ :, contrastNumber1, contrastNumber2 ]
            weights = np.sum( classSpecificWeights, 1 )
            
            # Build up stuff needed for rhs
            predicted = np.sum( classSpecificWeights * means[ :, contrastNumber2 ], 1 ) / ( weights + eps )
            residue = imageBuffers[ mask, contrastNumber2 ] - predicted
            tmp += weights * residue
            
            # Fill in submatrix of lhs
            weightsImageBuffer[ mask ] = weights
            lhs[ np.ix_( contrast1Indices, contrast2Indices ) ] \
                = computePrecisionOfKroneckerProductBasisFunctions( biasFieldBasisFunctions, 
                                                                    weightsImageBuffer )
            
        tmpImageBuffer[ mask ] = tmp
        rhs[ contrast1Indices ] = projectKroneckerProductBasisFunctions( biasFieldBasisFunctions, 
                                                                          tmpImageBuffer ).reshape(-1, 1)
          
    # Solve the linear system x = lhs \ rhs       
    solution = np.linalg.solve( lhs, rhs )
    
    #
    biasFieldCoefficients = solution.reshape( ( numberOfContrasts, numberOf3DBasisFunctions ) ).transpose()
    return biasFieldCoefficients



def deformMesh( mesh, transform, data, mask, means, variances, mixtureWeights, numberOfGaussiansPerClass, 
                userOptimizationParameters={} ):
  
    # Get images in ITK format
    numberOfContrasts = data.shape[ -1 ]
    images = []
    for contrastNumber in range( numberOfContrasts ):
        tmp = np.zeros( mask.shape, order='F' )
        tmp[ mask ] = data[ :, contrastNumber ]
        images.append( gems.KvlImage( requireNumpyArray( tmp ) ) )
        

    # Set up cost calculator
    calculator = gems.KvlCostAndGradientCalculator(
        typeName='AtlasMeshToIntensityImage',
        images=images,
        boundaryCondition='Sliding',
        transform=transform,
        means=means,
        variances=variances,
        mixtureWeights=mixtureWeights,
        numberOfGaussiansPerClass=numberOfGaussiansPerClass )

    # Get optimizer and plug calculator in it
    optimizerType = 'L-BFGS'
    optimizationParameters = {
        'Verbose': False,
        'MaximalDeformationStopCriterion': 0.001,  # measured in pixels,
        'LineSearchMaximalDeformationIntervalStopCriterion': 0.001,
        'MaximumNumberOfIterations': 20,
        'BFGS-MaximumMemoryLength': 12
         }
    optimizationParameters.update( userOptimizationParameters )
    print( optimizationParameters )
    optimizer = gems.KvlOptimizer( optimizerType, mesh, calculator, optimizationParameters )
    
    # Run deformation optimization
    historyOfDeformationCost = []
    historyOfMaximalDeformation = []
    nodePositionsBeforeDeformation = mesh.points
    while True:
        minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
        print( "maximalDeformation=%.4f minLogLikelihood=%.4f" % ( maximalDeformation, minLogLikelihoodTimesDeformationPrior ) )
        historyOfDeformationCost.append( minLogLikelihoodTimesDeformationPrior )
        historyOfMaximalDeformation.append( maximalDeformation )
        if maximalDeformation == 0:
            break
        
    # Return    
    nodePositionsAfterDeformation = mesh.points
    maximalDeformationApplied = np.sqrt(
        np.max( np.sum( ( nodePositionsAfterDeformation - nodePositionsBeforeDeformation) ** 2, 1 ) ) )
    return historyOfDeformationCost, historyOfMaximalDeformation, maximalDeformationApplied, minLogLikelihoodTimesDeformationPrior

  

  
def undoLogTransformAndBiasField( imageBuffers, biasFields, mask ):
    # 
    expBiasFields = np.zeros( biasFields.shape, order='F' )
    numberOfContrasts = imageBuffers.shape[-1]
    for contrastNumber in range( numberOfContrasts ):
        # We're computing it also outside of the mask, but clip the intensities there to the range 
        # observed inside the mask (with some margin) to avoid crazy extrapolation values
        biasField = biasFields[:, :, :, contrastNumber]
        clippingMargin = np.log( 2 )
        clippingMin = biasField[ mask ].min() - clippingMargin
        clippingMax = biasField[ mask ].max() + clippingMargin
        biasField[ biasField < clippingMin ] = clippingMin
        biasField[ biasField > clippingMax ] = clippingMax
        expBiasFields[ :, :, :, contrastNumber ] = np.exp( biasField )
        
    # 
    expImageBuffers = np.exp( imageBuffers ) / expBiasFields
    
    #
    return expImageBuffers, expBiasFields


def writeImage( fileName, buffer, cropping, example ):
  
  # Write un-cropped image to file
  uncroppedBuffer = np.zeros( example.getImageBuffer().shape, dtype=np.float32, order='F' )
  uncroppedBuffer[ cropping ] = buffer
  gems.KvlImage( requireNumpyArray( uncroppedBuffer ) ).write( fileName, example.transform_matrix )


def initializeGMMParameters( data, classPriors, numberOfGaussiansPerClass, useDiagonalCovarianceMatrices=True ):
  
    #
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    numberOfContrasts = data.shape[ -1 ]
    
  
    # Initialize the mixture parameters
    means = np.zeros( ( numberOfGaussians, numberOfContrasts ) )
    variances = np.zeros( ( numberOfGaussians, numberOfContrasts, numberOfContrasts ) )
    mixtureWeights = np.zeros( numberOfGaussians )
    for classNumber in range( numberOfClasses ):
        # Calculate the global weighted mean and variance of this class, where the weights are given by the prior
        prior = classPriors[ :, classNumber ]
        mean = data.T @ prior / np.sum( prior )
        tmp = data - mean
        prior = np.expand_dims( prior, 1 )
        variance = tmp.T @ ( tmp * prior ) / np.sum( prior )
        if useDiagonalCovarianceMatrices:
            # Force diagonal covariance matrices
            variance = np.diag( np.diag( variance ) )

        # Based on this, initialize the mean and variance of the individual Gaussian components in this class'
        # mixture model: variances are simply copied from the global class variance, whereas the means are
        # determined by splitting the [ mean-sqrt( variance ) mean+sqrt( variance ) ] domain into equal intervals,
        # the middle of which are taken to be the means of the Gaussians. Mixture weights are initialized to be
        # all equal.

        # This actually creates a mixture model that mimics the single Gaussian quite OK-ish
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]

        for componentNumber in range( numberOfComponents ):
            gaussianNumber = sum( numberOfGaussiansPerClass[ : classNumber ] ) + componentNumber
            variances[ gaussianNumber, :, : ] = variance
            intervalSize = 2 * np.sqrt( np.diag( variance ) ) / numberOfComponents
            means[ gaussianNumber, : ] = ( mean - np.sqrt( np.diag( variance ) ) + intervalSize / 2 + 
                                           componentNumber * intervalSize ).T
            mixtureWeights[ gaussianNumber ] = 1 / numberOfComponents

    #
    return means, variances, mixtureWeights

  
def logTransform( imageBuffers, mask ):

    logImageBuffers = imageBuffers.copy()
    logImageBuffers[ np.logical_not( mask ), : ] = 1
    logImageBuffers = np.log( logImageBuffers )

    #
    return logImageBuffers

  
def estimateModelParameters( imageBuffers, mask, biasFieldBasisFunctions, transform, voxelSpacing,
                             K, useDiagonalCovarianceMatrices,
                             classFractions, numberOfGaussiansPerClass, optimizationOptions,  
                             initialMeans=None, initialVariances=None, initialMixtureWeights=None,                              
                             initialBiasFieldCoefficients=None, initialDeformation=None, initialDeformationAtlasFileName=None,
                             hyperMeans=None, hyperMeansNumberOfMeasurements=None, 
                             hyperVariances=None, hyperVariancesNumberOfMeasurements=None,
                             hyperMixtureWeights=None, hyperMixtureWeightsNumberOfMeasurements=None, 
                             saveHistory=False, visualizer=None,
                             skipGMMParameterEstimationInFirstIteration=False,
                             skipBiasFieldParameterEstimationInFirstIteration=True,
                             hyperpriorPlugin=None
                             ):

    #  
    logger = logging.getLogger(__name__)
    history = []
    optimizationSummary = []
    if visualizer is None: visualizer = initVisualizer( False, False )
    
    # Convert optimizationOptions from dictionary into something more convenient to access
    optimizationOptions = Specification( optimizationOptions ) 
    source = optimizationOptions.multiResolutionSpecification
    optimizationOptions.multiResolutionSpecification = []
    for levelNumber in range( len( source ) ):
        optimizationOptions.multiResolutionSpecification.append( Specification( source[ levelNumber ] ) )
    print( '====================' )
    print( optimizationOptions )
    print( '====================' )


    # Parameter initialization. Deformation is encoded as node displacements (in template space) in a specific atlas
    means, variances, mixtureWeights = initialMeans, initialVariances, initialMixtureWeights
    biasFieldCoefficients = initialBiasFieldCoefficients
    deformation, deformationAtlasFileName = initialDeformation, initialDeformationAtlasFileName 


    # Loop over resolution levels
    fs.printPeakMemory( 'samsegment starting resolution loop' )
    numberOfMultiResolutionLevels = len( optimizationOptions.multiResolutionSpecification )
    for multiResolutionLevel in range( numberOfMultiResolutionLevels ):
        
        logger.debug( 'multiResolutionLevel=%d', multiResolutionLevel )
        visualizer.start_movie( window_id='Mesh deformation (level ' + str( multiResolutionLevel ) + ')', 
                                title='Mesh Deformation - the movie (level ' + str( multiResolutionLevel ) + ')' )
        
        
        maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification[multiResolutionLevel].maximumNumberOfIterations
        estimateBiasField = optimizationOptions.multiResolutionSpecification[multiResolutionLevel].estimateBiasField
        historyOfCost = [ 1/eps ]
        logger.debug( 'maximumNumberOfIterations: %d', maximumNumberOfIterations )
        
        # Downsample the images, the mask, the mesh, and the bias field basis functions (integer)
        logger.debug( 'Setting up downsampled model' )
        downSamplingFactors = np.uint32( np.round( optimizationOptions.multiResolutionSpecification[
                                                      multiResolutionLevel ].targetDownsampledVoxelSpacing / voxelSpacing ) )
        downSamplingFactors[ downSamplingFactors < 1 ] = 1
        downSampledImageBuffers, downSampledMask, downSampledMesh, downSampledInitialDeformationApplied, \
            downSampledTransform, downSampledBiasFieldBasisFunctions = \
            getDownSampledModel( imageBuffers, mask,
                                optimizationOptions.multiResolutionSpecification[multiResolutionLevel].atlasFileName,
                                K, transform, biasFieldBasisFunctions,
                                downSamplingFactors,
                                deformation, deformationAtlasFileName )
          
        # Also downsample the strength of the hyperprior, if any  
        downSampledHyperMeans, downSampledHyperMeansNumberOfMeasurements, \
            downSampledHyperVariances, downSampledHyperVariancesNumberOfMeasurements, \
            downSampledHyperMixtureWeights, downSampledHyperMixtureWeightsNumberOfMeasurements = \
            getDownsampledHyperparameters( hyperMeans, hyperMeansNumberOfMeasurements, 
                                           hyperVariances, hyperVariancesNumberOfMeasurements,
                                           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements,
                                           downSamplingFactors )
        
        
        # Save initial position at the start of this multi-resolution level
        initialNodePositions = downSampledMesh.points
        initialNodePositionsInTemplateSpace = mapPositionsFromSubjectToTemplateSpace( initialNodePositions, downSampledTransform )
          

        # Set priors in mesh to the merged (super-structure) ones
        mergedAlphas = gems.kvlMergeAlphas( downSampledMesh.alphas, classFractions )
        downSampledMesh.alphas = mergedAlphas

        
        #
        #visualizer.show( mesh=downSampledMesh, images=downSampledImageBuffers, window_id='samsegment em',
        #                 title='Mesh Registration (EM)' )
        visualizer.show( mesh=downSampledMesh, images=downSampledImageBuffers, 
                         window_id='Mesh deformation (level ' + str( multiResolutionLevel ) + ')',
                         title='Mesh Deformation (level ' + str( multiResolutionLevel ) + ')' )
    
    
        if saveHistory:
            levelHistory = { 'historyWithinEachIteration': [] }
        

        # Main iteration loop over both EM and deformation
        for iterationNumber in range( maximumNumberOfIterations ):
            
            logger.debug( 'iterationNumber=%d', iterationNumber )
            fs.printPeakMemory( 'samsegment resolution %d iteration %d' % (multiResolutionLevel, iterationNumber) )

            # Part I: estimate Gaussian mixture model parameters, as well as bias field parameters using EM.

            # Get the priors at the current mesh position
            tmp = downSampledMesh.rasterize_2( downSampledMask.shape, -1 )
            downSampledClassPriors = tmp[ downSampledMask ] / 65535

            # Initialize the model parameters if needed
            if means is None:
                means, variances, mixtureWeights = initializeGMMParameters( downSampledImageBuffers[ downSampledMask, : ], 
                                                                            downSampledClassPriors, numberOfGaussiansPerClass, 
                                                                            useDiagonalCovarianceMatrices )
                
                if hyperpriorPlugin is not None:
                    #
                    hyperMeans, hyperMeansNumberOfMeasurements, \
                        hyperVariances, hyperVariancesNumberOfMeasurements, \
                        hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements \
                        = hyperpriorPlugin( downSampledImageBuffers[ downSampledMask, : ],
                                            downSampledClassPriors, numberOfGaussiansPerClass, 
                                            voxelSpacing )
                    
                    downSampledHyperMeans, downSampledHyperMeansNumberOfMeasurements, \
                            downSampledHyperVariances, downSampledHyperVariancesNumberOfMeasurements, \
                            downSampledHyperMixtureWeights, downSampledHyperMixtureWeightsNumberOfMeasurements = \
                            getDownsampledHyperparameters( hyperMeans, hyperMeansNumberOfMeasurements, 
                                                           hyperVariances, hyperVariancesNumberOfMeasurements,
                                                           hyperMixtureWeights, hyperMixtureWeightsNumberOfMeasurements,
                                                           downSamplingFactors )
                
            if biasFieldCoefficients is None:
                numberOfBasisFunctions = [ functions.shape[1] for functions in downSampledBiasFieldBasisFunctions ]
                numberOfContrasts = downSampledImageBuffers.shape[-1]
                biasFieldCoefficients = np.zeros( ( np.prod( numberOfBasisFunctions ), numberOfContrasts ) )
            
            
            # Start EM iterations
            historyOfEMCost = [ 1/eps ]
            EMIterationNumber = 0
            while True:  
                logger.debug( 'EMIterationNumber=%d', EMIterationNumber )

                # Precompute intensities after bias field correction for later use (really only caching something that 
                # doesn't really figure in the model -- the real variable is biasFieldCoefficients)
                downSampledBiasFields = getBiasFields( biasFieldCoefficients, 
                                                      downSampledBiasFieldBasisFunctions,  downSampledMask )
                downSampledData = downSampledImageBuffers[ downSampledMask, : ] - downSampledBiasFields[ downSampledMask, : ]
                visualizer.show( image_list=[ downSampledBiasFields[ ...,i ] 
                                              for i in range( downSampledBiasFields.shape[-1] ) ], 
                                auto_scale=True, window_id='bias field', title='Bias Fields' )



                # E-step: compute the downSampledGaussianPosteriors based on the current parameters
                downSampledGaussianPosteriors, minLogLikelihood = getGaussianPosteriors( downSampledData, downSampledClassPriors, 
                                                                                         means, variances, mixtureWeights, 
                                                                                         numberOfGaussiansPerClass )
                    
                    
                # Compute the log-posterior of the model parameters, and check for convergence
                minLogGMMParametersPrior = evaluateMinLogPriorOfGMMParameters( means, variances, mixtureWeights, 
                                                                               numberOfGaussiansPerClass, 
                                                                               downSampledHyperMeans,
                                                                               downSampledHyperMeansNumberOfMeasurements,
                                                                               downSampledHyperVariances,
                                                                               downSampledHyperVariancesNumberOfMeasurements,
                                                                               downSampledHyperMixtureWeights,
                                                                               downSampledHyperMixtureWeightsNumberOfMeasurements )                 
                  
                historyOfEMCost.append( minLogLikelihood + minLogGMMParametersPrior )
                visualizer.plot( historyOfEMCost[ 1: ], window_id='history of EM cost', 
                                 title = 'History of EM Cost (level: ' + str( multiResolutionLevel ) 
                                         + ' iteration: ' + str( iterationNumber ) + ')' )
                EMIterationNumber += 1
                changeCostEMPerVoxel = ( historyOfEMCost[-2] - historyOfEMCost[-1] ) / downSampledData.shape[0]
                changeCostEMPerVoxelThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion
                if ( EMIterationNumber == 100 ) or ( changeCostEMPerVoxel < changeCostEMPerVoxelThreshold ):
                    # Converged
                    print('EM converged!')
                    break


                # M-step: update the model parameters based on the current posterior
                #
                # First the mixture model parameters
                if not( ( iterationNumber == 0 ) and skipGMMParameterEstimationInFirstIteration ):
                    means, variances, mixtureWeights = fitGMMParameters( downSampledData, downSampledGaussianPosteriors, 
                                                                        numberOfGaussiansPerClass,             
                                                                        useDiagonalCovarianceMatrices,
                                                                        downSampledHyperMeans,
                                                                        downSampledHyperMeansNumberOfMeasurements,
                                                                        downSampledHyperVariances,
                                                                        downSampledHyperVariancesNumberOfMeasurements,
                                                                        downSampledHyperMixtureWeights,
                                                                        downSampledHyperMixtureWeightsNumberOfMeasurements )
                 
                    
                # Now update the parameters of the bias field model.
                if ( estimateBiasField and 
                     not ( ( iterationNumber == 0 ) and skipBiasFieldParameterEstimationInFirstIteration ) ):
                    biasFieldCoefficients = fitBiasFieldParameters( downSampledImageBuffers, 
                                                                    downSampledGaussianPosteriors, means, variances,
                                                                    downSampledBiasFieldBasisFunctions, 
                                                                    downSampledMask )
                # End test if bias field update
                
            # End loop over EM iterations
            historyOfEMCost = historyOfEMCost[ 1: ]

            # Visualize the posteriors
            if hasattr( visualizer, 'show_flag' ):
                tmp = np.zeros( downSampledMask.shape + ( downSampledGaussianPosteriors.shape[-1], ) )
                tmp[ downSampledMask, : ] = downSampledGaussianPosteriors
                visualizer.show( probabilities=tmp, images=downSampledImageBuffers, window_id='EM Gaussian posteriors', 
                                 title= 'EM Gaussian posteriors (level: ' + str( multiResolutionLevel ) 
                                        + ' iteration: ' + str( iterationNumber ) + ')' )

            # Part II: update the position of the mesh nodes for the current mixture model and bias field parameter estimates
            optimizationParameters = {
                'Verbose': optimizationOptions.verbose,
                'MaximalDeformationStopCriterion': optimizationOptions.maximalDeformationStopCriterion,
                'LineSearchMaximalDeformationIntervalStopCriterion': optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion,
                'MaximumNumberOfIterations': optimizationOptions.maximumNumberOfDeformationIterations,
                'BFGS-MaximumMemoryLength': optimizationOptions.BFGSMaximumMemoryLength
            }
            historyOfDeformationCost, historyOfMaximalDeformation, maximalDeformationApplied, minLogLikelihoodTimesDeformationPrior = \
                deformMesh( downSampledMesh, downSampledTransform, downSampledData, downSampledMask, 
                            means, variances, mixtureWeights, numberOfGaussiansPerClass, optimizationParameters )

                

            # print summary of iteration
            print('iterationNumber: %d' % iterationNumber)
            print('maximalDeformationApplied: %.4f' % maximalDeformationApplied)
            print('=======================================================')
            #if hasattr( visualizer, 'show_flag' ):
            #    tmp = np.zeros( downSampledImageBuffers.shape ); tmp[ downSampledMask, : ] = downSampledData
            #    visualizer.show( mesh=downSampledMesh, images=tmp, window_id='samsegment mesh deformation (level ' + str( multiResolutionLevel ) + ')',
            #                     title='Mesh Registration (Deformation)' )
            visualizer.show( mesh=downSampledMesh, images=downSampledImageBuffers, 
                             window_id='Mesh deformation (level ' + str( multiResolutionLevel ) + ')',
                             title='Mesh Deformation (level ' + str( multiResolutionLevel ) + ')' )
    

            # Save history of the estimation
            if saveHistory:
                levelHistory['historyWithinEachIteration'].append( {
                    'historyOfEMCost': historyOfEMCost,
                    'mixtureWeights': mixtureWeights,
                    'means': means,
                    'variances': variances,
                    'biasFieldCoefficients': biasFieldCoefficients,
                    'historyOfDeformationCost': historyOfDeformationCost,
                    'historyOfMaximalDeformation': historyOfMaximalDeformation,
                    'maximalDeformationApplied': maximalDeformationApplied
                } )


            # Check for convergence
            historyOfCost.append( minLogLikelihoodTimesDeformationPrior + minLogGMMParametersPrior )
            visualizer.plot( historyOfCost[1:], window_id='history of cost (level ' + str( multiResolutionLevel ) + ')', 
                             title='History of Cost (level ' + str( multiResolutionLevel ) + ')' )
            previousCost = historyOfCost[-2]
            currentCost = historyOfCost[-1]
            costChange = previousCost - currentCost
            perVoxelDecrease = costChange / np.count_nonzero( downSampledMask )
            perVoxelDecreaseThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion
            if perVoxelDecrease < perVoxelDecreaseThreshold:
                break

        # End loop over coordinate descent optimization (intensity model parameters vs. atlas deformation)

        # Visualize the mesh deformation across iterations
        visualizer.show_movie( window_id='Mesh deformation (level ' + str( multiResolutionLevel ) + ')' )
        
        # Log the final per-voxel cost
        optimizationSummary.append( { 'numberOfIterations':  iterationNumber + 1,
                                      'perVoxelCost': currentCost / np.count_nonzero( downSampledMask ) } )

        # Get the final node positions
        finalNodePositions = downSampledMesh.points

        # Transform back in template space (i.e., undoing the affine registration that we applied), and save for later usage
        finalNodePositionsInTemplateSpace = mapPositionsFromSubjectToTemplateSpace( finalNodePositions, downSampledTransform )

        # Record deformation delta here in lieu of maintaining history
        deformation = finalNodePositionsInTemplateSpace - initialNodePositionsInTemplateSpace + downSampledInitialDeformationApplied
        deformationAtlasFileName = optimizationOptions.multiResolutionSpecification[ multiResolutionLevel ].atlasFileName

        # Save history of the estimation
        if saveHistory:
            levelHistory['downSamplingFactors'] = downSamplingFactors
            levelHistory['downSampledImageBuffers'] = downSampledImageBuffers
            levelHistory['downSampledMask'] = downSampledMask
            levelHistory['downSampledTransformMatrix'] = downSampledTransform.as_numpy_array
            levelHistory['initialNodePositions'] = initialNodePositions
            levelHistory['finalNodePositions'] = finalNodePositions
            levelHistory['initialNodePositionsInTemplateSpace'] = initialNodePositionsInTemplateSpace
            levelHistory['finalNodePositionsInTemplateSpace'] = finalNodePositionsInTemplateSpace
            levelHistory['historyOfCost'] = historyOfCost
            levelHistory['priorsAtEnd'] = downSampledClassPriors
            levelHistory['posteriorsAtEnd'] = downSampledGaussianPosteriors
            history.append( levelHistory )

    # End resolution level loop
  
  
    #
    return means, variances, mixtureWeights, biasFieldCoefficients, deformation, deformationAtlasFileName, optimizationSummary, history
           

def segment( imageBuffers, mask, transform, biasFieldBasisFunctions,
             atlasFileName, deformation, deformationAtlasFileName,
             means, variances, mixtureWeights, biasFieldCoefficients,
             numberOfGaussiansPerClass, classFractions,
             posteriorPlugin=None,
             posteriorPluginDictionary=None
           ):
    #
    fs.printPeakMemory( 'samsegment starting segmentation' )

    # Get the final mesh
    mesh = getMesh( atlasFileName, transform, 
                    initialDeformation=deformation, initialDeformationMeshCollectionFileName=deformationAtlasFileName )
      
    # Get the priors as dictated by the current mesh position
    priors = mesh.rasterize( imageBuffers.shape[ 0:3 ], -1 )
    priors = priors[ mask, : ]

    # Get bias field corrected data
    biasFields = getBiasFields( biasFieldCoefficients, biasFieldBasisFunctions )
    data = imageBuffers[ mask, : ] - biasFields[ mask, : ]

    # Compute the posterior distribution of the various structures
    if posteriorPlugin is None:
        posteriors = getPosteriors( data, priors, means, variances, mixtureWeights, numberOfGaussiansPerClass, classFractions )
    else:
        posteriors = posteriorPlugin( data, priors, means, variances, mixtureWeights, 
                                      numberOfGaussiansPerClass, classFractions,
                                      posteriorPluginDictionary )
        
    #
    estimatedNodePositions = mapPositionsFromSubjectToTemplateSpace( mesh.points, transform )

    #
    return posteriors, biasFields, estimatedNodePositions


def scaleBiasFields( biasFields, imageBuffers, mask, posteriors, targetIntensity=None, targetSearchStrings=None, names=None ):

    # Subtract a constant from the bias fields such that after bias field correction and exp-transform, the
    # average intensiy in the target structures will be targetIntensity
    if targetIntensity is not None:
        data = imageBuffers[ mask, : ] - biasFields[ mask, : ]
        targetWeights = np.zeros( data.shape[ 0 ] )
        for searchString in targetSearchStrings:
            for structureNumber, name in enumerate( names ):
                if searchString in name:
                    targetWeights += posteriors[ :, structureNumber ]
        offsets = np.log( targetIntensity ) - np.log( np.exp( data ).T @ targetWeights / np.sum( targetWeights ) )
        biasFields -= offsets.reshape( [ 1, 1, 1, biasFields.shape[-1] ] )
        
        # 
        scalingFactors = np.exp( offsets )
    else:
        scalingFactors = np.ones( imageBuffers.shape[-1] )
        
    return scalingFactors
  

def writeResults( imageFileNames, savePath, imageBuffers, mask, biasFields, posteriors, FreeSurferLabels, cropping,
                  targetIntensity=None, targetSearchStrings=None, names=None,
                  threshold=None, thresholdSearchString=None, savePosteriors=False
                  ):

    # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
    if threshold is not None:
        # Figure out the structure number of the special snowflake structure
        for structureNumber, name in enumerate( names ):
            if thresholdSearchString in name:
                thresholdStructureNumber = structureNumber
                break
    
        # Threshold
        print( 'thresholding posterior of ', names[ thresholdStructureNumber ], 'with threshold:', threshold )
        tmp = posteriors[ :, thresholdStructureNumber ].copy()
        posteriors[ :, thresholdStructureNumber ] = posteriors[ :, thresholdStructureNumber ] > threshold
        
        # Majority voting
        structureNumbers = np.array( np.argmax( posteriors, 1 ), dtype=np.uint32 )

        # Undo thresholding in posteriors
        posteriors[ :, thresholdStructureNumber ] = tmp
        
    else:
        # Majority voting
        structureNumbers = np.array( np.argmax( posteriors, 1 ), dtype=np.uint32 )
    
    
    
    freeSurferSegmentation = np.zeros( imageBuffers.shape[ 0:3 ], dtype=np.uint16 )
    FreeSurferLabels = np.array( FreeSurferLabels, dtype=np.uint16 )
    freeSurferSegmentation[ mask ] = FreeSurferLabels[ structureNumbers ]

    #
    scalingFactors = scaleBiasFields( biasFields, imageBuffers, mask, posteriors, targetIntensity, targetSearchStrings, names )

    # Get corrected intensities and bias field images in the non-log transformed domain
    expImageBuffers, expBiasFields = undoLogTransformAndBiasField( imageBuffers, biasFields, mask )


    # Write out various images - segmentation first
    exampleImage = gems.KvlImage( imageFileNames[ 0 ] )
    image_base_path, _ = os.path.splitext( imageFileNames[ 0 ] )
    _, scanName = os.path.split( image_base_path )
    writeImage( os.path.join( savePath, scanName + '_crispSegmentation.nii' ), freeSurferSegmentation, cropping, exampleImage )
    for contrastNumber, imageFileName in enumerate( imageFileNames ):
        image_base_path, _ = os.path.splitext( imageFileName )
        _, scanName = os.path.split( image_base_path )

        # Bias field
        writeImage( os.path.join( savePath, scanName + '_biasField.nii' ), expBiasFields[ ..., contrastNumber ], 
                    cropping, exampleImage )

        # Bias field corrected image
        writeImage( os.path.join( savePath, scanName + '_biasCorrected.nii' ), expImageBuffers[ ..., contrastNumber ], 
                    cropping, exampleImage )

        # Save a note indicating the scaling factor
        with open( os.path.join( savePath, scanName + '_scaling-factor.txt' ), 'w' ) as f:
            print( scalingFactors[ contrastNumber ], file=f )

    if savePosteriors:
        posteriorPath = os.path.join(savePath, 'posteriors')
        os.makedirs(posteriorPath, exist_ok=True)
        for i, name in enumerate(names):
            pvol = np.zeros(imageBuffers.shape[:3], dtype=np.float32)
            pvol[mask] = posteriors[:, i]
            writeImage(os.path.join(posteriorPath, name + '.nii'), pvol, cropping, exampleImage)

    # Compute volumes in mm^3
    volumeOfOneVoxel = np.abs( np.linalg.det( exampleImage.transform_matrix.as_numpy_array[ 0:3, 0:3 ] ) )
    volumesInCubicMm = ( np.sum( posteriors, axis=0 ) ) * volumeOfOneVoxel

    return volumesInCubicMm
  

def saveDeformedAtlas( originalAtlasFileName, deformedAtlasFileName, arg, applyAsDeformation=False ):

    #
    mesh_collection = gems.KvlMeshCollection()
    mesh_collection.read( originalAtlasFileName )
    if not applyAsDeformation:
        position = arg
    else:
        position = mesh_collection.reference_mesh.points + arg
    mesh_collection.reference_mesh.points = requireNumpyArray( position )
    mesh_collection.write( deformedAtlasFileName )
    
    #
    return



def samsegment( imageFileNames, atlasDir, savePath,
                transformedTemplateFileName=None, 
                userModelSpecifications={}, userOptimizationOptions={},
                visualizer=None, saveHistory=False, saveMesh=False,
                targetIntensity=None, targetSearchStrings=None,
                hyperpriorPlugin=None,
                posteriorPlugin=None,
                posteriorPluginVariables=None,
                threshold=None, thresholdSearchString=None, savePosteriors=False
                ):

    # Get full model specifications and optimization options (using default unless overridden by user) 
    modelSpecifications = getModelSpecifications( atlasDir, userModelSpecifications )
    optimizationOptions = getOptimizationOptions( atlasDir, userOptimizationOptions )
        

    # Print specifications
    print('##----------------------------------------------')
    print('              Samsegment Options')
    print('##----------------------------------------------')
    print('output directory:', savePath)
    print('input images:', ', '.join([imageFileName for imageFileName in imageFileNames]))
    print('transformed template:', transformedTemplateFileName)
    print('modelSpecifications:', modelSpecifications)
    print('optimizationOptions:', optimizationOptions)


    #
    modelSpecifications = Specification( modelSpecifications )

    # Setup a null visualizer if necessary
    if visualizer is None: visualizer = initVisualizer( False, False )

    # Make sure we can write in the target/results directory
    os.makedirs( savePath, exist_ok=True )


    # =======================================================================================
    #
    # Perform affine registration if needed
    #
    # =======================================================================================
    if transformedTemplateFileName is None:
        templateFileName = os.path.join( atlasDir, 'template.nii' )
        affineRegistrationMeshCollectionFileName = os.path.join( atlasDir, 'atlasForAffineRegistration.txt.gz' )
        _, transformedTemplateFileName, _ = registerAtlas( imageFileNames[ 0 ],
                                                           affineRegistrationMeshCollectionFileName,
                                                           templateFileName,
                                                           savePath,
                                                           visualizer=visualizer )




    # =======================================================================================
    #
    # Preprocessing (reading and masking of data)
    #
    # =======================================================================================

    # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
    # translation, rotation, scaling, and skewing) as well - this transformation will later be used
    # to initially transform the location of the atlas mesh's nodes into the coordinate system of the image.
    imageBuffers, transform, voxelSpacing, cropping = readCroppedImages( imageFileNames, transformedTemplateFileName )

    # Background masking: simply setting intensity values outside of a very rough brain mask to zero
    # ensures that they'll be skipped in all subsequent computations
    imageBuffers, mask = maskOutBackground( imageBuffers, modelSpecifications.atlasFileName, transform,
                                            modelSpecifications.brainMaskingSmoothingSigma, 
                                            modelSpecifications.brainMaskingThreshold )

    # Let's prepare for the bias field correction that is part of the imaging model. It assumes
    # an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
    # transform the data first.
    imageBuffers = logTransform( imageBuffers, mask )

    # Our bias model is a linear combination of a set of basis functions. We are using so-called "DCT-II" basis functions, 
    # i.e., the lowest few frequency components of the Discrete Cosine Transform.
    biasFieldBasisFunctions = getBiasFieldBasisFunctions( imageBuffers.shape[ 0:3 ], 
                                                          modelSpecifications.biasFieldSmoothingKernelSize / voxelSpacing )

    # Visualize some stuff
    if hasattr( visualizer, 'show_flag' ):
        visualizer.show( mesh=getMesh( modelSpecifications.atlasFileName, transform ), 
                        shape=imageBuffers.shape, 
                        window_id='samsegment mesh', title='Mesh', 
                        names=modelSpecifications.names, legend_width=350 )
        visualizer.show( images=imageBuffers, window_id='samsegment images', 
                        title='Samsegment Masked and Log-Transformed Contrasts' )

        import matplotlib.pyplot as plt   # avoid importing matplotlib by default
        plt.ion()
        f = plt.figure( 'Bias field basis functions' )
        for dimensionNumber in range(3):
            plt.subplot( 2, 2, dimensionNumber+1 )
            plt.plot( biasFieldBasisFunctions[ dimensionNumber ] )
        plt.draw()



    # =======================================================================================
    #
    # Parameter estimation
    #
    # =======================================================================================

    # The fact that we consider neuroanatomical structures as mixtures of "super"-structures for the purpose of model
    # parameter estimaton, but at the same time represent each of these super-structures with a mixture of Gaussians,
    # creates something of a messy situation when implementing this stuff. To avoid confusion, let's define a few
    # conventions that we'll closely follow in the code as follows:
    #
    #   - classNumber = 1 ... numberOfClasses  -> indexes a specific super-structure (there are numberOfClasses superstructures)
    #   - numberOfGaussiansPerClass            -> a numberOfClasses-dimensional vector that indicates the number of components
    #                                             in the Gaussian mixture model associated with each class
    #   - gaussianNumber = 1 .... numberOfGaussians  -> indexes a specific Gaussian distribution; there are
    #                                                   numberOfGaussians = sum( numberOfGaussiansPerClass ) of those in total
    #   - classFractions -> a numberOfClasses x numberOfStructures table indicating in each column the mixing weights of the
    #                       various classes in the corresponding structure
    numberOfGaussiansPerClass = [ param.numberOfComponents for param in modelSpecifications.sharedGMMParameters ]
    classFractions, _ = gems.kvlGetMergingFractionsTable( modelSpecifications.names, modelSpecifications.sharedGMMParameters )


    means, variances, mixtureWeights, biasFieldCoefficients, \
    deformation, deformationAtlasFileName, optimizationSummary, optimizationHistory = \
            estimateModelParameters( imageBuffers, mask, biasFieldBasisFunctions, transform, voxelSpacing,
                                    modelSpecifications.K, modelSpecifications.useDiagonalCovarianceMatrices,
                                    classFractions, numberOfGaussiansPerClass, optimizationOptions,
                                    saveHistory=saveHistory, visualizer=visualizer,
                                    hyperpriorPlugin=hyperpriorPlugin
                                    )


    # =======================================================================================
    #
    # Segment the data using the estimate model parameters, and write results out
    #
    # =======================================================================================

    # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
    # with all the original labels (instead of the reduced "super"-structure labels we created)
    posteriorPluginDictionary = {}
    if posteriorPluginVariables is not None:
        for variableName in posteriorPluginVariables: 
            posteriorPluginDictionary[ variableName ] = eval( variableName )
    posteriors, biasFields, nodePositions = segment( imageBuffers, mask, transform, biasFieldBasisFunctions,
                                                     modelSpecifications.atlasFileName, deformation, deformationAtlasFileName,
                                                     means, variances, mixtureWeights, biasFieldCoefficients,
                                                     numberOfGaussiansPerClass, classFractions,
                                                     posteriorPlugin,
                                                     posteriorPluginDictionary 
                                                    )
      
      

    # Write out segmentation and bias field corrected volumes
    volumesInCubicMm = writeResults( imageFileNames, savePath, imageBuffers, mask, biasFields,
                                     posteriors, modelSpecifications.FreeSurferLabels, cropping,
                                     targetIntensity, targetSearchStrings, modelSpecifications.names,
                                     threshold, thresholdSearchString, savePosteriors=savePosteriors
                                     )


    # Save the final mesh collection
    if saveMesh:
        print( 'Saving the final mesh in template space' )
        image_base_path, _ = os.path.splitext( imageFileNames[ 0 ] )
        _, scanName = os.path.split( image_base_path )
        deformedAtlasFileName = os.path.join( savePath, scanName + '_meshCollection.txt' )
        saveDeformedAtlas( modelSpecifications.atlasFileName , deformedAtlasFileName, nodePositions )


    # Save the history of the parameter estimation process
    if saveHistory:
        history = {'input': {
            'imageFileNames': imageFileNames,
            'transformedTemplateFileName': transformedTemplateFileName,
            'modelSpecifications': modelSpecifications,
            'optimizationOptions': optimizationOptions,
            'savePath': savePath
        }}
        history['imageBuffers'] = imageBuffers
        history['mask'] = mask
        history[ 'historyWithinEachMultiResolutionLevel' ] = optimizationHistory
        history[ "labels" ] = modelSpecifications.FreeSurferLabels
        history[ "names"] = modelSpecifications.names
        history[ "volumesInCubicMm"] = volumesInCubicMm
        history[ "optimizationSummary"] = optimizationSummary
        with open( os.path.join( savePath, 'history.p' ), 'wb' ) as file:
            pickle.dump( history, file, protocol=pickle.HIGHEST_PROTOCOL )


    return modelSpecifications.FreeSurferLabels, modelSpecifications.names, volumesInCubicMm, optimizationSummary



def generateSubjectSpecificTemplate( imageFileNamesList, savePath ):

    sstDir = os.path.join( savePath, 'sst' )
    os.makedirs( sstDir, exist_ok=True )

    sstFileNames = []
    for contrastNumber, contrastImageFileNames in enumerate( zip( *imageFileNamesList ) ):
      
        # Read in the various time point images, and compute the average
        numberOfTimepoints = len( contrastImageFileNames )
        image0 = gems.KvlImage( contrastImageFileNames[ 0 ] )
        imageBuffer = image0.getImageBuffer().copy()
        for timepointNumber in range( 1, numberOfTimepoints ):
            imageBuffer += gems.KvlImage( contrastImageFileNames[ timepointNumber ] ).getImageBuffer()
        imageBuffer /= numberOfTimepoints

        # Create an ITK image and write to disk
        sst = gems.KvlImage( requireNumpyArray( imageBuffer ) )
        sstFilename = os.path.join( sstDir, 'contrast' + str( contrastNumber ) + '_sst.mgz' )
        sst.write( sstFilename, image0.transform_matrix )

        #
        sstFileNames.append( sstFilename )

    #
    return sstFileNames
  
  


def samsegmentLongitudinal( imageFileNamesList, atlasDir, savePath,
                            userModelSpecifications={}, userOptimizationOptions={},
                            visualizer=None, saveHistory=False, 
                            targetIntensity=None, targetSearchStrings=None,
                            numberOfIterations=5,
                            strengthOfLatentGMMHyperprior=1.0,
                            strengthOfLatentDeformationHyperprior=10.0, 
                            saveSSTResults=True,
                            updateLatentMeans=True,
                            updateLatentVariances=True,
                            updateLatentMixtureWeights=True,
                            updateLatentDeformation=True,
                            initializeLatentDeformationToZero=False
                          ):

    """
    Longitudinal version of samsegment  

    The idea is based on the generative model in the paper 
    
      Iglesias, Juan Eugenio, et al. 
      Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using subject-specific atlases.
      Neuroimage 141 (2016): 542-555,

    in which a subject-specific atlas is obtained by generating a random warp from the usual population atlas, and
    subsequently each time point is again randomly warped from this subject-specific atlas. The intermediate 
    subject-specific atlas is effectively a latent variable in the model, and it's function is to encourage the 
    various time points to have atlas warps that are similar between themselves, without having to define a priori 
    what these warps should be similar to. In the implementation provided here, the stiffness of the first warp 
    (denoted by K0, as in the paper) is taken as the stiffness used in ordinary samseg, and the stiffness of the 
    second warp (denoted by K1) is controlled by the setting
    
      strengthOfLatentDeformationHyperprior
    
    so that K1 = strengthOfLatentDeformationHyperprior * K0. In the Iglesias paper, the setting 
    
      strengthOfLatentDeformationHyperprior = 1.0
      
    was used.  
    
    The overall idea is extended here by adding also latent variables encouraging corresponding Gaussian mixture 
    models across time points to be similar across time -- again without having to define a priori what exactly 
    they should look like. For given values of these latent variables, they effectively act has hyperparameters 
    on the mixture model parameters, the strength of which is controlled through the setting
    
      strengthOfLatentGMMHyperprior
      
    which weighs the relative strength of this hyperprior relative to the relevant (expected) data term in the 
    estimation procedures of the mixture model parameters at each time point. This aspect can be switched off by 
    setting 
    
      strengthOfLatentGMMHyperprior = 0.0
    
    NOTE: The general longitudinal pipeline in FreeSurfer 6.0 is described in the paper
    
      Reuter, Martin, et al. 
      Within-subject template estimation for unbiased longitudinal image analysis.
      Neuroimage 61.4 (2012): 1402-1418,
      
    which is based on the idea of retaining some temporal consistency by simply initializing the model fitting 
    procedure across all time points in exactly the same way. This is achieved by first creating a subject-specific 
    template that is subsequently analyzed and the result of which is then used as a (very good) initialization 
    in each time point. This behavior can be mimicked by setting
    
      initializeLatentDeformationToZero = True
      numberOfIterations = 1
      strengthOfLatentGMMHyperprior = 0.0
      strengthOfLatentDeformationHyperprior = 1.0
      
    in the implementation provided here.   
    """


    # Get full model specifications and optimization options (using default unless overridden by user) 
    modelSpecifications = getModelSpecifications( atlasDir, userModelSpecifications )
    optimizationOptions = getOptimizationOptions( atlasDir, userOptimizationOptions )
        
    #
    modelSpecifications = Specification( modelSpecifications )

    # Setup a null visualizer if necessary
    if visualizer is None: visualizer = initVisualizer( False, False )

    # Make sure we can write in the target/results directory
    os.makedirs( savePath, exist_ok=True )



  
    # =======================================================================================
    #
    # Construction and affine registration of subject-specific template (sst)
    #
    # =======================================================================================
  
    # Create the output folder
    os.makedirs( savePath, exist_ok=True )

    # Generate the subject specific template (sst)
    sstFileNames = generateSubjectSpecificTemplate( imageFileNamesList, savePath )

    # Affine atlas registration to sst
    templateFileName = os.path.join( atlasDir, 'template.nii' )
    affineRegistrationMeshCollectionFileName = os.path.join( atlasDir, 'atlasForAffineRegistration.txt.gz' )

    worldToWorldTransformMatrix, transformedTemplateFileName, _ \
        = registerAtlas( sstFileNames[ 0 ], affineRegistrationMeshCollectionFileName, templateFileName, savePath, 
                        visualizer=visualizer )



    # =======================================================================================
    #
    # Preprocessing (reading and masking of data)
    #
    # =======================================================================================

    sstImageBuffers, transform, voxelSpacing, cropping = readCroppedImages( sstFileNames, transformedTemplateFileName )

    imageBuffersList = []
    for imageFileNames in imageFileNamesList:
        imageBuffers, _, _, _ = readCroppedImages( imageFileNames, transformedTemplateFileName )
        imageBuffersList.append( imageBuffers )


    # Put everything in a big 4-D matrix to derive one consistent mask across all time points
    numberOfTimepoints = len( imageFileNamesList )
    imageSize = sstImageBuffers.shape[ :3 ]
    numberOfContrasts = sstImageBuffers.shape[ -1 ]
    combinedImageBuffers = np.zeros( imageSize + ( numberOfContrasts * ( 1 + numberOfTimepoints ), ) )
    combinedImageBuffers[ ..., 0:numberOfContrasts ] = sstImageBuffers
    for timepointNumber in range( numberOfTimepoints ):
        combinedImageBuffers[ ..., ( timepointNumber + 1 ) * numberOfContrasts : 
                                  ( timepointNumber + 2 ) * numberOfContrasts ] = imageBuffersList[ timepointNumber ]

    combinedImageBuffers, mask = maskOutBackground( combinedImageBuffers, modelSpecifications.atlasFileName, transform,
                                                    modelSpecifications.brainMaskingSmoothingSigma, 
                                                    modelSpecifications.brainMaskingThreshold )
    combinedImageBuffers = logTransform( combinedImageBuffers, mask )


    # Retrieve the masked sst and time points
    sstImageBuffers = combinedImageBuffers[ ..., 0:numberOfContrasts ]
    for timepointNumber in range( numberOfTimepoints ):
        imageBuffersList[ timepointNumber ] = combinedImageBuffers[ ..., ( timepointNumber + 1 ) * numberOfContrasts : 
                                                                        ( timepointNumber + 2 ) * numberOfContrasts ]

    visualizer.show( images=sstImageBuffers, title='sst' )
    for timepointNumber in range( numberOfTimepoints ):
        visualizer.show( images=imageBuffersList[ timepointNumber ], title='time point ' + str( timepointNumber ) )


    biasFieldBasisFunctions = getBiasFieldBasisFunctions( sstImageBuffers.shape[ 0:3 ], 
                                                          modelSpecifications.biasFieldSmoothingKernelSize / voxelSpacing )



    # =======================================================================================
    #
    # Parameter estimation for SST
    #
    # =======================================================================================
    numberOfGaussiansPerClass = [ param.numberOfComponents for param in modelSpecifications.sharedGMMParameters ]
    classFractions, _ = gems.kvlGetMergingFractionsTable( modelSpecifications.names, modelSpecifications.sharedGMMParameters )

    # 
    sstMeans, sstVariances, sstMixtureWeights, sstBiasFieldCoefficients, \
      sstDeformation, sstDeformationAtlasFileName, sstOptimizationSummary, sstOptimizationHistory = \
            estimateModelParameters( sstImageBuffers, mask, biasFieldBasisFunctions, transform, voxelSpacing,
                                    modelSpecifications.K, modelSpecifications.useDiagonalCovarianceMatrices,
                                    classFractions, numberOfGaussiansPerClass, optimizationOptions,  
                                    saveHistory=True, visualizer=visualizer )

    if hasattr( visualizer, 'show_flag' ):
        import matplotlib.pyplot as plt   # avoid importing matplotlib by default
        plt.ion()    
        sstBiasFields = getBiasFields( sstBiasFieldCoefficients, biasFieldBasisFunctions,  mask )
        sstData = sstImageBuffers[ mask, : ] - sstBiasFields[ mask, : ]
        axsList = []
        for contrastNumber in range( numberOfContrasts ):
            f = plt.figure()
            numberOfAxes = 2 + numberOfTimepoints
            numberOfRows = np.int( np.ceil( np.sqrt( numberOfAxes ) ) )
            numberOfColumns = np.int( np.ceil( numberOfAxes / numberOfRows ) )
            axs = f.subplots( numberOfRows, numberOfColumns, sharex=True )
            ax = axs.ravel()[ 0 ]
            _, bins, _ = ax.hist( sstImageBuffers[ mask, contrastNumber ], 100  )
            ax.grid()
            ax.set_title( 'sst before bias field correction' )
            ax = axs.ravel()[ 1 ]
            ax.hist( sstData[ :, contrastNumber ], bins )
            ax.grid()
            ax.set_title( 'sst after bias field correction' )
            for timepointNumber in range( numberOfTimepoints ):
                ax = axs.ravel()[ 2 + timepointNumber ]
                ax.hist( imageBuffersList[ timepointNumber ][ mask, contrastNumber ], bins )
                ax.grid()
                ax.set_title( 'time point ' + str( timepointNumber ) )
            axsList.append( axs )
        plt.draw()


    if saveHistory:
        history = { 
                  "sstMeans": sstMeans,
                  "sstVariances": sstVariances,
                  "sstMixtureWeights": sstMixtureWeights,
                  "sstBiasFieldCoefficients": sstBiasFieldCoefficients,
                  "sstDeformation": sstDeformation,
                  "sstDeformationAtlasFileName": sstDeformationAtlasFileName,
                  "sstOptimizationSummary": sstOptimizationSummary,
                  "sstOptimizationHistory": sstOptimizationHistory
                  }


    if saveSSTResults:
        sstPosteriors, sstBiasFields, _ = segment( sstImageBuffers, 
                                                  mask, transform, biasFieldBasisFunctions,
                                                  modelSpecifications.atlasFileName, 
                                                  sstDeformation, 
                                                  sstDeformationAtlasFileName,
                                                  sstMeans, 
                                                  sstVariances, 
                                                  sstMixtureWeights,
                                                  sstBiasFieldCoefficients, 
                                                  numberOfGaussiansPerClass, classFractions )

        #
        sstDir, _ = os.path.split( sstFileNames[ 0 ] )

        #
        sstVolumesInCubicMm = writeResults( sstFileNames, 
                                            sstDir, sstImageBuffers, 
                                            mask, sstBiasFields, sstPosteriors, 
                                            modelSpecifications.FreeSurferLabels, cropping,
                                            targetIntensity, targetSearchStrings, modelSpecifications.names )

        #
        if saveHistory:
            history[ "sstVolumesInCubicMm" ] = sstVolumesInCubicMm




    # =======================================================================================
    #
    # Iterative parameter vs. latent variables estimation, using SST result for initialization
    # and/or anchoring of hyperprior strength
    #
    # =======================================================================================

    # Initialization of the time-specific model parameters
    timepointMeans, timepointVariances, timepointMixtureWeights, timepointBiasFieldCoefficients, \
        timepointDeformations, timepointDeformationAtlasFileNames = \
        [ sstMeans ] * numberOfTimepoints, [ sstVariances ] * numberOfTimepoints, [ sstMixtureWeights ] * numberOfTimepoints, \
        [ None ] * numberOfTimepoints, [ None ] * numberOfTimepoints, [ None ] * numberOfTimepoints



    # Initialization of the latent variables, acting as hyperparameters when viewed from the model parameters' perspective
    latentDeformation = sstDeformation.copy()
    latentDeformationAtlasFileName = sstDeformationAtlasFileName
    latentMeans = sstMeans.copy()
    latentVariances = sstVariances.copy()
    latentMixtureWeights = sstMixtureWeights.copy()

    if initializeLatentDeformationToZero:
        timepointDeformations, timepointDeformationAtlasFileNames = \
            [ latentDeformation ] * numberOfTimepoints, [ latentDeformationAtlasFileName ] * numberOfTimepoints
        latentDeformation[ : ] = 0



    # Strength of the hyperprior (i.e., how much the latent variables controll the conditional posterior of the parameters) 
    # is user-controlled. 
    # 
    # For the GMM part, I'm using the *average* number of voxels assigned to the components in each mixture (class) of the 
    # SST segmentation, so that all the components in each mixture are well-regualized (and tiny components don't get to do
    # whatever they want)
    K0 = modelSpecifications.K # Stiffness population -> latent position
    K1 = strengthOfLatentDeformationHyperprior * K0 # Stiffness latent position -> each time point
    sstEstimatedNumberOfVoxelsPerGaussian = np.sum( sstOptimizationHistory[ -1 ][ 'posteriorsAtEnd' ], axis=0 ) * \
                                              np.prod( sstOptimizationHistory[ -1 ][ 'downSamplingFactors' ] )
    numberOfClasses = len( numberOfGaussiansPerClass )
    numberOfGaussians = sum( numberOfGaussiansPerClass )
    latentMeansNumberOfMeasurements = np.zeros( numberOfGaussians ) 
    latentVariancesNumberOfMeasurements = np.zeros( numberOfGaussians )
    latentMixtureWeightsNumberOfMeasurements = np.zeros( numberOfClasses )
    for classNumber in range( numberOfClasses ):
        #
        numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
        gaussianNumbers = np.array( np.sum( numberOfGaussiansPerClass[ :classNumber ] ) + \
                                    np.array( range( numberOfComponents ) ), dtype=np.uint32 )
        sstEstimatedNumberOfVoxelsInClass = \
            np.sum( sstEstimatedNumberOfVoxelsPerGaussian[ gaussianNumbers ] )

        latentMixtureWeightsNumberOfMeasurements[ classNumber ] = strengthOfLatentGMMHyperprior * sstEstimatedNumberOfVoxelsInClass
      
        averageSizeOfComponents = sstEstimatedNumberOfVoxelsInClass / numberOfComponents
        latentMeansNumberOfMeasurements[ gaussianNumbers ] = strengthOfLatentGMMHyperprior * averageSizeOfComponents
        latentVariancesNumberOfMeasurements[ gaussianNumbers ] = strengthOfLatentGMMHyperprior * averageSizeOfComponents


    # Estimating the mode of the latentVariance posterior distribution (which is Wishart) requires a stringent condition 
    # on latentVariancesNumberOfMeasurements so that the mode is actually defined
    threshold = ( numberOfContrasts + 1 ) + 1e-6
    latentVariancesNumberOfMeasurements[ latentVariancesNumberOfMeasurements < threshold ] = threshold

    # No point in updating latent GMM parameters if the GMM hyperprior has zero weight. The latent variances are also
    # a bit tricky, as they're technically driven to zero in that scenario -- let's try not to go there...
    if ( strengthOfLatentGMMHyperprior == 0 ):
        updateLatentMeans, updateLatentVariances, updateLatentMixtureWeights = False, False, False


    # Loop over all iterations
    historyOfTotalCost, historyOfTotalTimepointCost, historyOfLatentAtlasCost = [], [], []
    progressPlot = None
    iterationNumber = 0
    if saveHistory:
        history = { **history, 
                    **{
                      "timepointMeansEvolution": [],
                      "timepointVariancesEvolution": [],
                      "timepointMixtureWeightsEvolution": [],
                      "timepointBiasFieldCoefficientsEvolution": [],
                      "timepointDeformationsEvolution": [],
                      "timepointDeformationAtlasFileNamesEvolution": [],
                      "latentMeansEvolution": [],
                      "latentVariancesEvolution": [],
                      "latentMixtureWeightsEvolution": [],
                      "latentDeformationEvolution": [],
                      "latentDeformationAtlasFileNameEvolution": []
                      }
                  } 
        
    while True:

        # =======================================================================================
        #
        # Update parameters for each time point using the current latent variable estimates
        #
        # =======================================================================================

        # Create a new atlas that will be the basis to deform the individual time points from
        latentAtlasFileName = os.path.join( savePath, 'latentAtlas_iteration' + str( iterationNumber ) )
        saveDeformedAtlas( latentDeformationAtlasFileName, latentAtlasFileName, latentDeformation, True )

        # Only use the last resolution level, and with the newly created atlas as atlas
        timepointOptimizationOptions = optimizationOptions.copy()
        timepointOptimizationOptions[ 'multiResolutionSpecification' ] = [ timepointOptimizationOptions[ 'multiResolutionSpecification' ][-1] ]
        timepointOptimizationOptions[ 'multiResolutionSpecification' ][ 0 ][ 'atlasFileName' ] = latentAtlasFileName
        print( timepointOptimizationOptions )

        # Loop over all time points
        totalTimepointCost = 0
        for timepointNumber in range( numberOfTimepoints ):
            # 
            timepointMeans[ timepointNumber ], timepointVariances[ timepointNumber ], timepointMixtureWeights[ timepointNumber ], \
                timepointBiasFieldCoefficients[ timepointNumber ], \
                timepointDeformations[ timepointNumber ], timepointDeformationAtlasFileNames[ timepointNumber ], \
                optimizationSummary, optimizationHistory = \
                    estimateModelParameters( imageBuffersList[ timepointNumber ], mask, biasFieldBasisFunctions, transform, voxelSpacing,
                                              K1, modelSpecifications.useDiagonalCovarianceMatrices,
                                              classFractions, numberOfGaussiansPerClass, timepointOptimizationOptions,  
                                              saveHistory=True, visualizer=visualizer,
                                              initialMeans=timepointMeans[ timepointNumber ], 
                                              initialVariances=timepointVariances[ timepointNumber ], 
                                              initialMixtureWeights=timepointMixtureWeights[ timepointNumber ],
                                              initialBiasFieldCoefficients=timepointBiasFieldCoefficients[ timepointNumber ],
                                              initialDeformation=timepointDeformations[ timepointNumber ], 
                                              initialDeformationAtlasFileName=timepointDeformationAtlasFileNames[ timepointNumber ],
                                              hyperMeans=latentMeans, 
                                              hyperMeansNumberOfMeasurements=latentMeansNumberOfMeasurements, 
                                              hyperVariances=latentVariances,
                                              hyperVariancesNumberOfMeasurements=latentVariancesNumberOfMeasurements,
                                              hyperMixtureWeights=latentMixtureWeights, 
                                              hyperMixtureWeightsNumberOfMeasurements=latentMixtureWeightsNumberOfMeasurements,
                                              skipBiasFieldParameterEstimationInFirstIteration=False,
                                              skipGMMParameterEstimationInFirstIteration=( iterationNumber == 0 )
                                            )
              
            totalTimepointCost += optimizationHistory[ -1 ][ 'historyOfCost' ][ -1 ]

            print( '=================================' )
            print( '\n' )
            print( 'timepointNumber: ', timepointNumber )
            print( 'perVoxelCost: ', optimizationSummary[-1][ 'perVoxelCost' ] )
            print( '\n' )
            print( '=================================' )
            if hasattr( visualizer, 'show_flag' ):
                import matplotlib.pyplot as plt   # avoid importing matplotlib by default
                plt.ion()    
                timepointBiasFields = getBiasFields( timepointBiasFieldCoefficients[ timepointNumber ], 
                                                     biasFieldBasisFunctions,  mask )
                timepointData = imageBuffersList[ timepointNumber ][ mask, : ] - timepointBiasFields[ mask, : ]
                for contrastNumber in range( numberOfContrasts ):
                    axs = axsList[ contrastNumber ]
                    ax = axs.ravel()[ 2 + timepointNumber ]
                    ax.clear()
                    ax.hist( timepointData[ :, contrastNumber ], bins )
                    ax.grid()
                    ax.set_title( 'time point ' + str( timepointNumber ) )
                plt.draw()

            # End loop over time points
            
            


        # =======================================================================================
        #
        # Check for convergence. 
        # =======================================================================================

        # In order to also measure the deformation from the population atlas -> latent position,
        # create:
        #   (1) a mesh collection with as reference position the population reference position, and as positions 
        #       the currently estimated time point positions. 
        #   (2) a mesh with the current latent position 
        # Note that in (1) we don't need those time positions now, but these will come in handy very soon to 
        # optimize the latent position
        # 
        # The parameter estimation happens in a (potentially) downsampled image grid, so it's import to work in the same space 
        # when measuring and updating the latentDeformation
        transformUsedForEstimation = gems.KvlTransform( requireNumpyArray( sstOptimizationHistory[ -1 ][ 'downSampledTransformMatrix' ] ) )
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read( latentDeformationAtlasFileName )
        mesh_collection.transform( transformUsedForEstimation )
        referencePosition = mesh_collection.reference_position
        timepointPositions = []
        for timepointNumber in range( numberOfTimepoints ):
            positionInTemplateSpace = mapPositionsFromSubjectToTemplateSpace( referencePosition, transformUsedForEstimation ) + \
                                        latentDeformation + timepointDeformations[ timepointNumber ]
            timepointPositions.append( mapPositionsFromTemplateToSubjectSpace( positionInTemplateSpace, transformUsedForEstimation ) )
        mesh_collection.set_positions( referencePosition, timepointPositions )

        # Read mesh in sst warp
        mesh = getMesh( latentAtlasFileName, transformUsedForEstimation )
        
        # 
        calculator = gems.KvlCostAndGradientCalculator( mesh_collection, K0, 0.0, transformUsedForEstimation )
        latentAtlasCost, _ = calculator.evaluate_mesh_position( mesh )
    
        #
        totalCost = totalTimepointCost + latentAtlasCost
        print( '*' * 100 + '\n' )
        print( 'iterationNumber: ', iterationNumber )
        print( 'totalCost: ', totalCost )
        print( '   latentAtlasCost: ', latentAtlasCost )
        print( '   totalTimepointCost: ', totalTimepointCost )
        print( '*' * 100 + '\n' )
        historyOfTotalCost.append( totalCost ) 
        historyOfTotalTimepointCost.append( totalTimepointCost )
        historyOfLatentAtlasCost.append( latentAtlasCost )
        
        if hasattr( visualizer, 'show_flag' ):
            import matplotlib.pyplot as plt   # avoid importing matplotlib by default
            plt.ion()
            if progressPlot is None:
                plt.figure()
                progressPlot = plt.subplot()
            progressPlot.clear()
            progressPlot.plot( historyOfTotalCost, color='k' )
            progressPlot.plot( historyOfTotalTimepointCost, linestyle='-.', color='b' )
            progressPlot.plot( historyOfLatentAtlasCost, linestyle='-.', color='r'  )
            progressPlot.grid()
            progressPlot.legend( [ 'total', 'timepoints', 'latent atlas deformation' ] )
            plt.draw()

        if saveHistory:
            history[ "timepointMeansEvolution" ].append( timepointMeans.copy() )
            history[ "timepointVariancesEvolution" ].append( timepointVariances.copy() )
            history[ "timepointMixtureWeightsEvolution" ].append( timepointMixtureWeights.copy() )
            history[ "timepointBiasFieldCoefficientsEvolution" ].append( timepointBiasFieldCoefficients.copy() )
            history[ "timepointDeformationsEvolution" ].append( timepointDeformations.copy() )
            history[ "timepointDeformationAtlasFileNamesEvolution" ].append( timepointDeformationAtlasFileNames.copy() )
            history[ "latentMeansEvolution" ].append( latentMeans.copy() )
            history[ "latentVariancesEvolution" ].append( latentVariances.copy() )
            history[ "latentMixtureWeightsEvolution" ].append( latentMixtureWeights.copy() )
            history[ "latentDeformationEvolution" ].append( latentDeformation.copy() )
            history[ "latentDeformationAtlasFileNameEvolution" ].append( latentDeformationAtlasFileName )

        if iterationNumber >= ( numberOfIterations-1 ):
            print( 'Stopping' )
            break

            

        # =======================================================================================
        #
        # Update the latent variables based on the current parameter estimates
        #
        # =======================================================================================

        # Update the latentDeformation
        if updateLatentDeformation:
            # Set up calculator
            calculator = gems.KvlCostAndGradientCalculator( mesh_collection, K0, K1, transformUsedForEstimation )

            # Get optimizer and plug calculator in it
            optimizerType = 'L-BFGS'
            optimizationParameters = {
                'Verbose': optimizationOptions[ 'verbose' ],
                'MaximalDeformationStopCriterion': optimizationOptions[ 'maximalDeformationStopCriterion' ],
                'LineSearchMaximalDeformationIntervalStopCriterion': optimizationOptions[ 'lineSearchMaximalDeformationIntervalStopCriterion' ],
                'MaximumNumberOfIterations': optimizationOptions[ 'maximumNumberOfDeformationIterations' ],
                'BFGS-MaximumMemoryLength': optimizationOptions[ 'BFGSMaximumMemoryLength' ]
            }
            optimizer = gems.KvlOptimizer( optimizerType, mesh, calculator, optimizationParameters )



            # Run deformation optimization
            historyOfDeformationCost = []
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            while True:
                minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
                print( "maximalDeformation=%.4f minLogLikelihood=%.4f" % ( maximalDeformation, minLogLikelihoodTimesDeformationPrior ) )
                historyOfDeformationCost.append( minLogLikelihoodTimesDeformationPrior )
                historyOfMaximalDeformation.append( maximalDeformation )
                if maximalDeformation == 0:
                    break
                
            nodePositionsAfterDeformation = mesh.points
            maximalDeformationApplied = np.sqrt(
                np.max( np.sum( ( nodePositionsAfterDeformation - nodePositionsBeforeDeformation) ** 2, 1 ) ) )


            # 
            nodePositionsBeforeDeformation = mapPositionsFromSubjectToTemplateSpace( nodePositionsBeforeDeformation, transformUsedForEstimation )
            nodePositionsAfterDeformation = mapPositionsFromSubjectToTemplateSpace( nodePositionsAfterDeformation, transformUsedForEstimation )
            estimatedUpdate = nodePositionsAfterDeformation - nodePositionsBeforeDeformation
            latentDeformation += estimatedUpdate


        # Update latentMeans
        if updateLatentMeans:
            numberOfGaussians =  np.sum( numberOfGaussiansPerClass )
            numberOfContrasts = latentMeans.shape[ -1 ]
            for gaussianNumber in range( numberOfGaussians ):
                # Set up linear system
                lhs = np.zeros( ( numberOfContrasts, numberOfContrasts ) )
                rhs = np.zeros( ( numberOfContrasts, 1 ) )
                for timepointNumber in range( numberOfTimepoints ):
                    mean = np.expand_dims( timepointMeans[ timepointNumber ][ gaussianNumber ], 1 )
                    variance = timepointVariances[ timepointNumber ][ gaussianNumber ]
                    
                    lhs += np.linalg.inv( variance )
                    rhs += np.linalg.solve( variance, mean )
                
                # Solve linear system
                latentMean = np.linalg.solve( lhs, rhs )
                latentMeans[ gaussianNumber, : ] = latentMean.T


        # Update latentVariances
        if updateLatentVariances:
            numberOfGaussians =  np.sum( numberOfGaussiansPerClass )
            numberOfContrasts = latentMeans.shape[ -1 ]
            for gaussianNumber in range( numberOfGaussians ):
                # Precision is essentially averaged 
                averagePrecision = np.zeros( ( numberOfContrasts, numberOfContrasts ) )
                for timepointNumber in range( numberOfTimepoints ):
                    variance = timepointVariances[ timepointNumber ][ gaussianNumber ]
                    
                    averagePrecision += np.linalg.inv( variance )
                averagePrecision /= numberOfTimepoints
                
                latentVarianceNumberOfMeasurements = latentVariancesNumberOfMeasurements[ gaussianNumber ]
                latentVariance = np.linalg.inv( averagePrecision ) * \
                                ( latentVarianceNumberOfMeasurements - numberOfContrasts - 1 ) / latentVarianceNumberOfMeasurements
                latentVariances[ gaussianNumber ] = latentVariance

          
        # Update latentMixtureWeights
        if updateLatentMixtureWeights:
            numberOfClasses = len( numberOfGaussiansPerClass )
            for classNumber in range( numberOfClasses ):
                numberOfComponents = numberOfGaussiansPerClass[ classNumber ]
                averageInLogDomain = np.zeros( numberOfComponents )
                for componentNumber in range( numberOfComponents ):
                    gaussianNumber = sum( numberOfGaussiansPerClass[ :classNumber ] ) + componentNumber
                    for timepointNumber in range( numberOfTimepoints ):
                        mixtureWeight = timepointMixtureWeights[ timepointNumber ][ gaussianNumber ]
                        averageInLogDomain[ componentNumber ] += np.log( mixtureWeight + eps )
                    averageInLogDomain[ componentNumber ] /= numberOfTimepoints
                
                # Solution is normalized version
                solution = np.exp( averageInLogDomain )
                solution /= np.sum( solution + eps )
                
                #
                for componentNumber in range( numberOfComponents ):
                    gaussianNumber = sum( numberOfGaussiansPerClass[ :classNumber ] ) + componentNumber
                    latentMixtureWeights[ gaussianNumber ] = solution[ componentNumber ]

          
        iterationNumber += 1
        # End loop over parameter and latent variable estimation iterations



    # =======================================================================================
    #
    # Using estimated parameters, segment and write out results for each time point
    #
    # =======================================================================================
    timepointVolumesInCubicMm = []
    for timepointNumber in range( numberOfTimepoints ): 
        #
        posteriors, biasFields, nodePositions = segment( imageBuffersList[ timepointNumber ], 
                                                        mask, transform, biasFieldBasisFunctions,
                                                        modelSpecifications.atlasFileName, 
                                                        latentDeformation + timepointDeformations[ timepointNumber ], 
                                                        latentDeformationAtlasFileName,
                                                        timepointMeans[ timepointNumber ], 
                                                        timepointVariances[ timepointNumber ], 
                                                        timepointMixtureWeights[ timepointNumber ],
                                                        timepointBiasFieldCoefficients[ timepointNumber ], 
                                                        numberOfGaussiansPerClass, classFractions )

        #
        timepointDir = os.path.join( savePath, 'timepoint' + str( timepointNumber ) )
        os.makedirs( timepointDir, exist_ok=True )
        
        #
        volumesInCubicMm = writeResults( imageFileNamesList[ timepointNumber ], 
                                        timepointDir, imageBuffersList[ timepointNumber ], 
                                        mask, biasFields, posteriors, 
                                        modelSpecifications.FreeSurferLabels, cropping,
                                        targetIntensity, targetSearchStrings, modelSpecifications.names )

        timepointVolumesInCubicMm.append( volumesInCubicMm )



    #
    optimizationSummary = { "historyOfTotalCost": historyOfTotalCost, 
                            "historyOfTotalTimepointCost": historyOfTotalTimepointCost, 
                            "historyOfLatentAtlasCost": historyOfLatentAtlasCost }

    if saveHistory:
        history[ "labels" ] = modelSpecifications.FreeSurferLabels
        history[ "names"] = modelSpecifications.names
        history[ "timepointVolumesInCubicMm"] = timepointVolumesInCubicMm
        history[ "optimizationSummary"] = optimizationSummary
        with open( os.path.join( savePath, 'historyLongitudinal.p' ), 'wb' ) as file:
            pickle.dump( history, file, protocol=pickle.HIGHEST_PROTOCOL )

    return modelSpecifications.FreeSurferLabels, modelSpecifications.names, timepointVolumesInCubicMm, optimizationSummary



