import os
import numpy as np
import itertools
import freesurfer as fs

from .utilities import icv
from .io import kvlReadCompressionLookupTable, kvlReadSharedGMMParameters
from .figures import initVisualizer
from .utilities import requireNumpyArray
from . import gemsbindings as gems


def getModelSpecifications(atlasDir, userModelSpecifications={}, pallidumAsWM=True, gmmFileName=None):

    # Create default model specifications as a dictionary
    FreeSurferLabels, names, colors = kvlReadCompressionLookupTable(os.path.join(atlasDir, 'compressionLookupTable.txt'))

    # Use default sharedGMMParameters.txt file in the atlas directory if no custom file is provided
    if gmmFileName is None:
        gmmFileName = os.path.join(atlasDir, 'sharedGMMParameters.txt')
    else:
        # If the custom GMM file does not exist in the working dir, assume it exists in the atlas dir
        if not os.path.isfile(gmmFileName):
            gmmFileName = os.path.join(atlasDir, gmmFileName)
    if not os.path.isfile(gmmFileName):
        fs.fatal('GMM parameter file does not exist at %s' % gmmFileName)

    sharedGMMParameters = kvlReadSharedGMMParameters(gmmFileName)

    # If pallidumAsWM is True remove from the sharedGMMParameters 'Pallidum' as an independent class
    # and move it into 'GlobalWM'.
    if pallidumAsWM:
        pallidumGMMNumber = None
        globalWMGMMNumber = None
        for classNumber, mergeOption in enumerate(sharedGMMParameters):
            if 'Pallidum' == mergeOption.mergedName:
                pallidumGMMNumber = classNumber
            elif 'GlobalWM' == mergeOption.mergedName:
                globalWMGMMNumber = classNumber

        if pallidumGMMNumber is not None and globalWMGMMNumber is not None:
            sharedGMMParameters[globalWMGMMNumber].searchStrings.append('Pallidum')
            sharedGMMParameters.pop(pallidumGMMNumber)

    modelSpecifications = {
        'FreeSurferLabels': FreeSurferLabels,
        'atlasFileName': os.path.join(atlasDir, 'atlas_level2.txt.gz'),
        'names': names,
        'colors': colors,
        'sharedGMMParameters': sharedGMMParameters,
        'useDiagonalCovarianceMatrices': True,
        'maskingProbabilityThreshold': 0.5, # threshold on probability of background
        'maskingDistance': 10.0, # distance in mm of how far into background the mask goes out
        'K': 0.1,  # stiffness of the mesh
        'biasFieldSmoothingKernelSize': 50,  # distance in mm of sinc function center to first zero crossing
    }

    modelSpecifications.update(userModelSpecifications)

    return modelSpecifications


def getOptimizationOptions(atlasDir, userOptimizationOptions={}):

    # Create default optimization options as a dictionary
    optimizationOptions = {
        'maximumNumberOfDeformationIterations': 20,
        'absoluteCostPerVoxelDecreaseStopCriterion': 1e-4,
        'verbose': False,
        'maximalDeformationStopCriterion': 0.001,  # measured in pixels
        'lineSearchMaximalDeformationIntervalStopCriterion': 0.001,
        'maximalDeformationAppliedStopCriterion': 0.0,
        'BFGSMaximumMemoryLength': 12,
        'multiResolutionSpecification': [
            {
                # level 1
                'atlasFileName': os.path.join(atlasDir, 'atlas_level1.txt.gz'),
                'targetDownsampledVoxelSpacing': 2.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True
            }, {
                # level 2
                'atlasFileName': os.path.join(atlasDir, 'atlas_level2.txt.gz'),
                'targetDownsampledVoxelSpacing': 1.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True
            }
        ]
    }

    # Overwrite with any user specified options. The 'multiResolutionSpecification' key has as value a list
    # of dictionaries which we shouldn't just over-write, but rather update themselves, so this is special case
    userOptimizationOptionsCopy = userOptimizationOptions.copy()
    key = 'multiResolutionSpecification'
    if key in userOptimizationOptionsCopy:
        userList = userOptimizationOptionsCopy[key]
        defaultList = optimizationOptions[key]
        for levelNumber in range(len(defaultList)):
            if levelNumber < len(userList):
                defaultList[levelNumber].update(userList[levelNumber])
            else:
                del defaultList[levelNumber]
        del userOptimizationOptionsCopy[key]
    optimizationOptions.update(userOptimizationOptionsCopy)

    return optimizationOptions


def readCroppedImages(imageFileNames, templateFileName, imageToImageTransform):
    # Read the image data from disk and crop it given a template image and it's associated
    # registration matrix.

    croppedImageBuffers = []
    for imageFileName in imageFileNames:

        input_image = fs.Volume.read(imageFileName)
        template_image = fs.Volume.read(templateFileName)

        imageToImage = fs.LinearTransform(imageToImageTransform)

        # Map each of the corners of the bounding box, and record minima and maxima
        boundingLimit = np.array(template_image.shape[:3]) - 1
        corners = np.array(list(itertools.product(*zip((0, 0, 0), boundingLimit))))
        transformedCorners = imageToImage.transform(corners)

        inputLimit = np.array(input_image.shape[:3]) - 1
        minCoord = np.clip(transformedCorners.min(axis=0).astype(int),     (0, 0, 0), inputLimit)
        maxCoord = np.clip(transformedCorners.max(axis=0).astype(int) + 1, (0, 0, 0), inputLimit) + 1

        cropping = tuple([slice(min, max) for min, max in zip(minCoord, maxCoord)])
        croppedImageBuffers.append(input_image.data[cropping])

        # create and translate kvl transform
        transform = imageToImage.matrix.copy()
        transform[:3, 3] -= minCoord
        transform = gems.KvlTransform(requireNumpyArray(transform))

    croppedImageBuffers = np.transpose(croppedImageBuffers, axes=[1, 2, 3, 0])

    # Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels,
    # downsampling steps etc in mm.
    nonCroppedImage = gems.KvlImage(imageFileNames[0])
    imageToWorldTransformMatrix = nonCroppedImage.transform_matrix.as_numpy_array
    voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)

    return croppedImageBuffers, transform, voxelSpacing, cropping


def readCroppedImagesLegacy(imageFileNames, transformedTemplateFileName):
    # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
    # translation, rotation, scaling, and skewing) as well - this transformation will later be used
    # to initially transform the location of the atlas mesh's nodes into the coordinate system of the image.
    imageBuffers = []
    for imageFileName in imageFileNames:
        # Get the pointers to image and the corresponding transform
        image = gems.KvlImage(imageFileName, transformedTemplateFileName)
        transform = image.transform_matrix
        cropping = image.crop_slices
        imageBuffers.append(image.getImageBuffer())

    imageBuffers = np.transpose(imageBuffers, axes=[1, 2, 3, 0])

    # Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels,
    # downsampling steps etc in mm.
    nonCroppedImage = gems.KvlImage(imageFileNames[0])
    imageToWorldTransformMatrix = nonCroppedImage.transform_matrix.as_numpy_array
    voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)

    #
    return imageBuffers, transform, voxelSpacing, cropping


def showImage(data):
    range = (data.min(), data.max())

    Nx = data.shape[0]
    Ny = data.shape[1]
    Nz = data.shape[2]

    x = round(Nx / 2)
    y = round(Ny / 2)
    z = round(Nz / 2)

    xySlice = data[:, :, z]
    xzSlice = data[:, y, :]
    yzSlice = data[x, :, :]

    patchedSlices = np.block([[xySlice, xzSlice], [yzSlice.T, np.zeros((Nz, Nz)) + range[0]]])

    import matplotlib.pyplot as plt  # avoid importing matplotlib by default
    plt.imshow(patchedSlices.T, cmap=plt.cm.gray, vmin=range[0], vmax=range[1])
    # plt.gray()
    # plt.imshow( patchedSlices.T, vmin=range[ 0 ], vmax=range[ 1 ] )
    # plt.show()
    plt.axis('off')


def maskOutBackground(imageBuffers, atlasFileName, transform, 
                      maskingProbabilityThreshold, maskingDistance,
                      probabilisticAtlas, voxelSpacing, visualizer=None, maskOutZeroIntensities=True):
    # Setup a null visualizer if necessary
    if visualizer is None:
        visualizer = initVisualizer(False, False)

    # Read the affinely coregistered atlas mesh (in reference position)
    mesh = probabilisticAtlas.getMesh(atlasFileName, transform)

    # Mask away uninteresting voxels. This is done by a poor man's implementation of a dilation operation on
    # a non-background class mask; followed by a cropping to the area covered by the mesh (needed because
    # otherwise there will be voxels in the data with prior probability zero of belonging to any class)
    imageSize = imageBuffers.shape[0:3]
    labelNumber = 0
    backgroundPrior = mesh.rasterize_1a(imageSize, labelNumber)

    if os.environ.get('SAMSEG_LEGACY_BACKGROUND_MASKING') is not None:
        print('INFO: using legacy background masking option')
        
        #
        brainMaskingSmoothingSigma = 3.0
        brainMaskingThreshold = 0.01
        
        # Threshold background prior at 0.5 - this helps for atlases built from imperfect (i.e., automatic)
        # segmentations, whereas background areas don't have zero probability for non-background structures
        backGroundThreshold = 2 ** 8
        backGroundPeak = 2 ** 16 - 1
        backgroundPrior = np.ma.filled(np.ma.masked_greater(backgroundPrior, backGroundThreshold),
                                      backGroundPeak).astype(np.float32)

        visualizer.show(probabilities=backgroundPrior, images=imageBuffers, window_id='samsegment background',
                        title='Background Priors')

        smoothingSigmas = [1.0 * brainMaskingSmoothingSigma] * 3
        smoothedBackgroundPrior = gems.KvlImage.smooth_image_buffer(backgroundPrior, smoothingSigmas)
        visualizer.show(probabilities=smoothedBackgroundPrior, window_id='samsegment smoothed',
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
        brainMaskThreshold = 65535.0 * (1.0 - brainMaskingThreshold)
        mask = np.ma.less(smoothedBackgroundPrior, brainMaskThreshold)
        
    else:
        #
        #visualizer = initVisualizer(True, True)
        from scipy import ndimage
      
        # Threshold prior of background
        backgroundMask = np.ma.greater( backgroundPrior, (2**16-1) * maskingProbabilityThreshold )
        visualizer.show( images=backgroundMask.astype(float), title='thresholded' )

        # Extend by distance of maskingDistance (in mm)
        distance = ndimage.distance_transform_edt( backgroundMask, sampling=voxelSpacing )
        print( voxelSpacing )
        mask = np.ma.less( distance, maskingDistance )
        visualizer.show( images=mask.astype(float), title='short distance' )

        # Fill holes inside the mask, if any        
        mask = ndimage.binary_fill_holes( mask ) 
        visualizer.show( images=mask.astype(float), title='holes filled' )
        

    # Crop to area covered by the mesh
    alphas = mesh.alphas
    areaCoveredAlphas = [[0.0, 1.0]] * alphas.shape[0]
    mesh.alphas = areaCoveredAlphas  # temporary replacement of alphas
    areaCoveredByMesh = mesh.rasterize_1b(imageSize, 1)
    mesh.alphas = alphas  # restore alphas
    mask = np.logical_and(mask, areaCoveredByMesh)

    # If a pixel has a zero intensity in any of the contrasts, that is also masked out across all contrasts
    if maskOutZeroIntensities:
        numberOfContrasts = imageBuffers.shape[-1]
        for contrastNumber in range(numberOfContrasts):
            mask *= imageBuffers[:, :, :, contrastNumber] > 0

    # Mask the images
    maskedImageBuffers = imageBuffers.copy()
    maskedImageBuffers[np.logical_not(mask), :] = 0

    #
    return maskedImageBuffers, mask


def undoLogTransformAndBiasField(imageBuffers, biasFields, mask):
    #
    expBiasFields = np.zeros(biasFields.shape, order='F')
    numberOfContrasts = imageBuffers.shape[-1]
    for contrastNumber in range(numberOfContrasts):
        # We're computing it also outside of the mask, but clip the intensities there to the range
        # observed inside the mask (with some margin) to avoid crazy extrapolation values
        biasField = biasFields[:, :, :, contrastNumber]
        clippingMargin = np.log(2)
        clippingMin = biasField[mask].min() - clippingMargin
        clippingMax = biasField[mask].max() + clippingMargin
        biasField[biasField < clippingMin] = clippingMin
        biasField[biasField > clippingMax] = clippingMax
        expBiasFields[:, :, :, contrastNumber] = np.exp(biasField)

    #
    expImageBuffers = np.exp(imageBuffers) / expBiasFields
    expImageBuffers[ np.logical_not( mask ), : ] = 0
    
    #
    return expImageBuffers, expBiasFields


def writeImage(fileName, buffer, cropping, example):

    # Write un-cropped image to file
    uncroppedBuffer = np.zeros(example.getImageBuffer().shape, dtype=np.float32, order='F')
    uncroppedBuffer[cropping] = buffer
    gems.KvlImage(requireNumpyArray(uncroppedBuffer)).write(fileName, example.transform_matrix)


def logTransform(imageBuffers, mask):

    logImageBuffers = imageBuffers.copy().astype( 'float' )
    logImageBuffers[ logImageBuffers == 1 ] += 1e-5 # Voxels with zero values but inside the mask 
                                                    # should not be skipped in the C++ code!
    logImageBuffers[np.logical_not(mask), :] = 1
    logImageBuffers = np.log(logImageBuffers)

    #
    return logImageBuffers


def scaleBiasFields(biasFields, imageBuffers, mask, posteriors, targetIntensity=None, targetSearchStrings=None,
                        names=None):

    # Subtract a constant from the bias fields such that after bias field correction and exp-transform, the
    # average intensiy in the target structures will be targetIntensity
    if targetIntensity is not None:
        data = imageBuffers[mask, :] - biasFields[mask, :]
        targetWeights = np.zeros(data.shape[0])
        for searchString in targetSearchStrings:
            for structureNumber, name in enumerate(names):
                if searchString in name:
                    targetWeights += posteriors[:, structureNumber]
        offsets = np.log(targetIntensity) - np.log(np.exp(data).T @ targetWeights / np.sum(targetWeights))
        biasFields -= offsets.reshape([1, 1, 1, biasFields.shape[-1]])

        #
        scalingFactors = np.exp(offsets)
    else:
        scalingFactors = np.ones(imageBuffers.shape[-1])

    return scalingFactors


def convertRASTransformToLPS(ras2ras):
    ras2lps = np.diag([-1, -1, 1, 1])
    return ras2lps @ ras2ras @ np.linalg.inv(ras2lps)


def convertLPSTransformToRAS(lps2lps):
    ras2lps = np.diag([-1, -1, 1, 1])
    return np.linalg.inv(ras2lps) @ lps2lps @ ras2lps
