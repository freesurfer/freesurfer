import os
import scipy.io
import math
import numpy as np
import logging
import pickle
from functools import reduce
from operator import mul

import freesurfer as fs
import freesurfer.gems as gems

from .utilities import Specification, requireNumpyArray, ensureDims
from .figures import initVisualizer
from .bias_correction import projectKroneckerProductBasisFunctions, backprojectKroneckerProductBasisFunctions, \
    computePrecisionOfKroneckerProductBasisFunctions, biasCorrectData


logger = logging.getLogger(__name__)
eps = np.finfo(float).eps


def samsegment(
    imageFileNames,
    transformedTemplateFileName,
    modelSpecifications,
    optimizationOptions,
    savePath,
    visualizer=None,
    saveHistory=False,
    saveMesh=False,
):

    # Print specifications
    print('##----------------------------------------------')
    print('              Samsegment Options')
    print('##----------------------------------------------')
    print('output directory:', savePath)
    print('input images:', ', '.join([imageFileName for imageFileName in imageFileNames]))
    print('transformed template:', transformedTemplateFileName)
    print('modelSpecifications:', modelSpecifications)
    print('optimizationOptions:', optimizationOptions)

    # Save input variables in a history dictionary
    if saveHistory:
        history = {'input': {
            'imageFileNames': imageFileNames,
            'transformedTemplateFileName': transformedTemplateFileName,
            'modelSpecifications': modelSpecifications,
            'optimizationOptions': optimizationOptions,
            'savePath': savePath
        }}

    # Setup a null visualizer if necessary
    if visualizer is None: visualizer = initVisualizer(False, False)

    # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
    # translation, rotation, scaling, and skewing) as well - this transformation will later be used
    # to initially transform the location of the atlas mesh's nodes into the coordinate system of the image.

    imageBuffers = []
    images = []
    for imageFileName in imageFileNames:
        # Get the pointers to image and the corresponding transform
        image = gems.KvlImage(imageFileName, transformedTemplateFileName)
        transform = image.transform_matrix
        nonCroppedImageSize = image.non_cropped_image_size
        cropping = image.crop_slices
        images.append(image)
        imageBuffers.append(image.getImageBuffer())

    imageSize = imageBuffers[0].shape
    imageBuffers = np.transpose(imageBuffers, axes=[1, 2, 3, 0])
    visualizer.show(images=imageBuffers, window_id='samsegment contrast', title='Samsegment Contrasts')

    # Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels, downsampling steps etc in mm.
    nonCroppedImage = gems.KvlImage(imageFileNames[0])
    imageToWorldTransformMatrix = nonCroppedImage.transform_matrix.as_numpy_array
    voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)

    # Read the atlas mesh from file, immediately applying the previously determined transform to the location
    # of its nodes. Rather than one mesh, the atlas consists of a so-called "collection" of meshes: they're
    # all the same except for the location of their mesh nodes. The first mesh, so-called "reference mesh"
    # has index -1: it represents the average shape (location of the mesh nodes) of all the other meshes. If
    # the atlas was built from N training subjects, there will be N other meshes as well, each one warping
    # the atlas optimally to each of the training subjects

    # You also need to provide a value for K, which determines the flexibility of the atlas mesh, i.e., how
    # much it will typically deform. Higher values correspond to stiffer meshes.
    meshCollection = gems.KvlMeshCollection()
    meshCollection.read(modelSpecifications.atlasFileName)
    meshCollection.k = modelSpecifications.K
    meshCollection.transform(transform)

    # Retrieve the reference mesh, i.e., the mesh representing the average shape.
    mesh = meshCollection.reference_mesh

    # Get a Matlab matrix containing a copy of the probability vectors in each mesh node (size numberOfNodes x
    # numberOfLabels ).
    alphas = mesh.alphas

    # Mask away uninteresting voxels. This is done by a poor man's implementation of a dilation operation on
    # a non-background class mask; followed by a cropping to the area covered by the mesh (needed because
    # otherwise there will be voxels in the data with prior probability zero of belonging to any class)
    labelNumber = 0
    backgroundPrior = mesh.rasterize_1a(imageSize, labelNumber)

    # Threshold background prior at 0.5 - this helps for atlases built from imperfect (i.e., automatic)
    # segmentations, whereas background areas don't have zero probability for non-background structures
    backGroundThreshold = 2 ** 8
    backGroundPeak = 2 ** 16 - 1
    backgroundPrior = np.ma.filled(
        np.ma.masked_greater(backgroundPrior, backGroundThreshold),
        backGroundPeak).astype(np.float32)

    visualizer.show(probabilities=backgroundPrior, images=imageBuffers, window_id='samsegment background', title='Background Priors')
    smoothingSigmas = [1.0 * modelSpecifications.brainMaskingSmoothingSigma] * 3
    smoothedBackgroundPrior = gems.KvlImage.smooth_image_buffer(backgroundPrior, smoothingSigmas)
    visualizer.show(probabilities=smoothedBackgroundPrior, window_id='samsegment smoothed', title='Smoothed Background Priors')

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
    brainMaskThreshold = 65535.0 * (1.0 - modelSpecifications.brainMaskingThreshold)
    brainMask = np.ma.less(smoothedBackgroundPrior, brainMaskThreshold)

    # Crop to area covered by the mesh
    areaCoveredAlphas = [[0.0, 1.0]] * alphas.shape[0]
    mesh.alphas = areaCoveredAlphas  # temporary replacement of alphas
    areaCoveredByMesh = mesh.rasterize_1b(imageSize, 1)
    mesh.alphas = alphas  # restore alphas
    brainMask = np.logical_and(brainMask, areaCoveredByMesh)

    # Mask each of the inputs
    numberOfContrasts = len(imageFileNames)
    for contrastNumber in range(numberOfContrasts):
        imageBuffers[:, :, :, contrastNumber] *= brainMask

    visualizer.show(images=imageBuffers, window_id='samsegment images', title='Samsegment Masked Contrasts')

    # Let's prepare for the bias field correction that is part of the imaging model. It assumes
    # an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
    # transform the data first. In order to do so, mask out zeros from the images.
    # This removes any voxel where any contrast has a zero value (messes up log)
    mask = np.full(imageSize, True, dtype=np.bool)
    for contrastNumber in range(numberOfContrasts):
        mask *= imageBuffers[:, :, :, contrastNumber] > 0
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings('ignore')
        log_buffers = np.log(imageBuffers)

    imageBuffers = np.ma.fix_invalid(log_buffers).filled(0)
    for contrastNumber in range(numberOfContrasts):
        imageBuffers[np.logical_not(mask), contrastNumber] = 0
    log_buffers = None

    if saveHistory:
        history['imageBuffers'] = imageBuffers
        history['mask'] = mask

    # Merge classes into "super-structures" that define a single Gaussian mixture model shared between the classes belonging
    # to the same super-structure
    FreeSurferLabels = modelSpecifications.FreeSurferLabels
    names = modelSpecifications.names
    colors = modelSpecifications.colors
    [reducedAlphas, reducedNames, reducedFreeSurferLabels, reducedColors, translationTable
        ] = gems.kvlMergeAlphas(alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors)

    visualizer.show(mesh=mesh, shape=imageBuffers.shape, window_id='samsegment mesh', title='Mesh', names=names, legend_width=350)

    # The fact that we merge several neuroanatomical structures into "super"-structures for the purpose of model
    # parameter estimaton, but at the same time represent each of these super-structures with a mixture of Gaussians,
    # creates something of a messy situation when implementing this stuff. To avoid confusion, let's define a few
    # conventions that we'll closely follow in the code as follows:
    #
    #   - classNumber = 1 ... numberOfClasses  -> indexes a specific super-structure (there are numberOfClasses superstructures)
    #   - numberOfGaussiansPerClass            -> a numberOfClasses-dimensional vector that indicates the number of components
    #                                             in the Gaussian mixture model associated with each class
    #   - gaussianNumber = 1 .... numberOfGaussians  -> indexes a specific Gaussian distribution; there are
    #                                                   numberOfGaussians = sum( numberOfGaussiansPerClass ) of those in total

    # In certain situations it will be easier to index things using "classNumber" (which is the only relevant thing when
    # estimating mesh deformations) than using "gaussianNumber" (which is the only relevant thing when estimating mixture
    # model parameters) and vice versa, so we need to have a standardized way of going back and forth between those. For a
    # given classNumber, there are numberOfComponents = numberOfGaussiansPerClass( classNumber ) components in its mixture
    # model. By convention we convert the pair ( classNumber, componentNumber ) into gaussianNumber as follows:
    numberOfGaussiansPerClass = [param.numberOfComponents for param in modelSpecifications.sharedGMMParameters]
    numberOfClasses = len(numberOfGaussiansPerClass)
    numberOfGaussians = sum(numberOfGaussiansPerClass)

    # Our bias model is a linear combination of a set of basis functions. We are using so-called
    # "DCT-II" basis functions, i.e., the lowest few frequency components of the Discrete Cosine
    # Transform.
    kroneckerProductBasisFunctions = []
    numberOfBasisFunctions = []
    for dimensionNumber in range(3):
        N = imageSize[dimensionNumber]
        delta = modelSpecifications.biasFieldSmoothingKernelSize / voxelSpacing[dimensionNumber]
        M = math.ceil(N / delta) + 1
        Nvirtual = (M - 1) * delta
        js = [(index + 0.5) * math.pi / Nvirtual for index in range(N)]
        scaling = [math.sqrt(2 / Nvirtual)] * M
        scaling[0] /= math.sqrt(2)
        A = np.array([[math.cos(freq * m) * scaling[m] for m in range(M)] for freq in js])
        kroneckerProductBasisFunctions.append(A)
        numberOfBasisFunctions.append(M)
    basisProduct = reduce(mul, numberOfBasisFunctions, 1)
    biasFieldCoefficients = np.zeros((basisProduct, numberOfContrasts))

    if saveHistory:
        history['historyWithinEachMultiResolutionLevel'] = []

    fs.printPeakMemory('samsegment starting resolution loop')

    numberOfMultiResolutionLevels = len(optimizationOptions.multiResolutionSpecification)
    for multiResolutionLevel in range(numberOfMultiResolutionLevels):
        logger.debug('multiResolutionLevel=%d', multiResolutionLevel)
        #  If the movie flag is on then making a movie archives a lot of data.
        #  Saving some memory here by making, showing, then erasing the movie at each resolution level.
        visualizer.start_movie(window_id='samsegment', title='Samsegment Mesh Registration - the movie')
        maximumNumberOfIterations = optimizationOptions.multiResolutionSpecification[multiResolutionLevel].maximumNumberOfIterations
        estimateBiasField = optimizationOptions.multiResolutionSpecification[multiResolutionLevel].estimateBiasField
        historyOfCost = [1 / eps]
        logger.debug('maximumNumberOfIterations: %d', maximumNumberOfIterations)
        # Downsample the images, the mask, the mesh, and the bias field basis functions
        # Must be integer
        downSamplingFactors = np.uint32(np.round(optimizationOptions.multiResolutionSpecification[
                                                     multiResolutionLevel].targetDownsampledVoxelSpacing / voxelSpacing))
        downSamplingFactors[downSamplingFactors < 1] = 1
        downSampledMask = mask[::downSamplingFactors[0], ::downSamplingFactors[1], ::downSamplingFactors[2]]
        downSampledMaskIndices = np.where(downSampledMask)
        activeVoxelCount = len(downSampledMaskIndices[0])
        downSampledImageBuffers = np.zeros(downSampledMask.shape + (numberOfContrasts,), order='F')
        for contrastNumber in range(numberOfContrasts):
            logger.debug('first time contrastNumber=%d', contrastNumber)
            # No image smoothing
            downSampledImageBuffers[:, :, :, contrastNumber] = imageBuffers[::downSamplingFactors[0],
                                                               ::downSamplingFactors[1],
                                                               ::downSamplingFactors[2],
                                                               contrastNumber]

        downSampledKroneckerProductBasisFunctions = [np.array(kroneckerProductBasisFunction[::downSamplingFactor])
                                                     for kroneckerProductBasisFunction, downSamplingFactor in
                                                     zip(kroneckerProductBasisFunctions, downSamplingFactors)]
        downSampledImageSize = downSampledImageBuffers[:, :, :, 0].shape
        # Read the atlas mesh to be used for this multi-resolution level, taking into account
        # the downsampling to position it correctly
        downSamplingTransformMatrix = np.diag(1. / downSamplingFactors)
        downSamplingTransformMatrix = np.pad(downSamplingTransformMatrix, (0, 1), mode='constant', constant_values=0)
        downSamplingTransformMatrix[3][3] = 1

        totalTransformationMatrix = downSamplingTransformMatrix @ transform.as_numpy_array

        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(optimizationOptions.multiResolutionSpecification[multiResolutionLevel].atlasFileName)
        mesh_collection.k = modelSpecifications.K
        mesh_collection.transform(gems.KvlTransform(requireNumpyArray(totalTransformationMatrix)))

        mesh = mesh_collection.reference_mesh

        # Get the initial mesh node positions, also transforming them back into template space
        # (i.e., undoing the affine registration that we applied) for later usage
        initialNodePositions = mesh.points
        numberOfNodes = len(initialNodePositions)
        tmp = np.linalg.solve(totalTransformationMatrix,
                              np.pad(initialNodePositions, ((0, 0), (0, 1)), mode='constant', constant_values=1).T).T
        initialNodePositionsInTemplateSpace = tmp[:, 0:3]
        # If this is not the first multi-resolution level, apply the warp computed during the previous level
        if multiResolutionLevel > 0:
            # Get the warp in template space
            [initialNodeDeformationInTemplateSpace, initial_averageDistance, initial_maximumDistance] = gems.kvlWarpMesh(
                optimizationOptions.multiResolutionSpecification[multiResolutionLevel - 1].atlasFileName,
                nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel,
                optimizationOptions.multiResolutionSpecification[multiResolutionLevel].atlasFileName)
            # Apply this warp on the mesh node positions in template space, and transform into current space
            desiredNodePositionsInTemplateSpace = initialNodePositionsInTemplateSpace + initialNodeDeformationInTemplateSpace
            tmp = (totalTransformationMatrix @ np.pad(desiredNodePositionsInTemplateSpace, ((0, 0), (0, 1)), 'constant', constant_values=1).T).T
            desiredNodePositions = tmp[:, 0:3]
            mesh.points = requireNumpyArray(desiredNodePositions)

        # Set priors in mesh to the reduced (super-structure) ones
        alphas = mesh.alphas
        reducedAlphas, _, _, _, _ = gems.kvlMergeAlphas(alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors)
        mesh.alphas = reducedAlphas

        # Algorithm-wise, we're just estimating sets of parameters for one given data (MR scan) that is
        # known and fixed throughout. However, in terms of bias field correction it will be computationally
        # more efficient to pre-compute the bias field corrected version of the scan ("corrected" with
        # the current estimate of the bias field) once and pass that on to different routines instead of the
        # original data.
        # For convenience (although potentially a recipe for future bug introduction), I'm also keeping a
        # vectorized form of that around -- this will be useful in various places in the EM-parts. So
        # effectively I have two redundant variables "downSampledBiasCorrectedImageBuffers" and "biasCorrectedData"
        # that really just encode the variable "biasFieldCoefficients" and so need to be meticiously updated each time
        # "biasFieldCoefficients" is updated (!)
        downSampledBiasCorrectedImageBuffers = np.zeros(downSampledImageSize + (numberOfContrasts,), order='F')
        biasCorrectedData = np.zeros((activeVoxelCount, numberOfContrasts), order='F')

        downSampledBiasFields = biasCorrectData(
            biasCorrectedData,
            biasFieldCoefficients,
            downSampledBiasCorrectedImageBuffers,
            downSampledImageBuffers,
            downSampledKroneckerProductBasisFunctions,
            downSampledMask,
            downSampledMaskIndices,
            numberOfContrasts
        )
        visualizer.show( image_list=downSampledBiasFields, auto_scale=True, window_id='bias field', title='Bias Fields')

        # Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
        # we start deforming. We'll use this just for visualization purposes
        posteriors = np.zeros((activeVoxelCount, numberOfGaussians), order='F')

        # Easier to work with vector notation in the EM computations
        # reshape into a matrix
        data = np.zeros((activeVoxelCount, numberOfContrasts))
        for contrastNumber in range(numberOfContrasts):
            tmp = downSampledImageBuffers[:, :, :, contrastNumber]
            data[:, contrastNumber] = tmp[downSampledMaskIndices]

        if saveHistory:
            levelHistory = {'historyWithinEachIteration': []}

        # Main iteration loop over both EM and deformation
        for iterationNumber in range(maximumNumberOfIterations):
            logger.debug('iterationNumber=%d', iterationNumber)
            fs.printPeakMemory('samsegment resolution %d iteration %d' % (multiResolutionLevel, iterationNumber))

            # Part I: estimate Gaussian mixture model parameters, as well as bias field parameters using EM.

            # Get the priors at the current mesh position
            tmp = mesh.rasterize_2(downSampledImageSize, -1)
            priors = tmp[downSampledMaskIndices] / 65535

            # Start EM iterations
            if ((multiResolutionLevel == 0) and (iterationNumber == 0)):
                #
                # Initialize the mixture parameters if this is the first time ever you run this
                means = np.zeros((numberOfGaussians, numberOfContrasts))
                variances = np.zeros((numberOfGaussians, numberOfContrasts, numberOfContrasts))
                mixtureWeights = np.zeros((numberOfGaussians, 1))
                for classNumber in range(numberOfClasses):
                    # Calculate the global weighted mean and variance of this class, where the weights are given by the prior
                    prior = priors[:, classNumber]
                    mean = data.T @ prior / np.sum(prior)
                    tmp = data - mean
                    prior = np.expand_dims(prior, 1)
                    variance = tmp.T @ (tmp * prior) / np.sum(prior)
                    if modelSpecifications.useDiagonalCovarianceMatrices:
                        # Force diagonal covariance matrices
                        variance = np.diag(np.diag(variance))

                    # Based on this, initialize the mean and variance of the individual Gaussian components in this class'
                    # mixture model: variances are simply copied from the global class variance, whereas the means are
                    # determined by splitting the [ mean-sqrt( variance ) mean+sqrt( variance ) ] domain into equal intervals,
                    # the middle of which are taken to be the means of the Gaussians. Mixture weights are initialized to be
                    # all equal.

                    # This actually creates a mixture model that mimics the single Gaussian quite OK-ish: to visualize this do e.g.
                    numberOfComponents = numberOfGaussiansPerClass[classNumber]

                    for componentNumber in range(numberOfComponents):
                        gaussianNumber = sum(numberOfGaussiansPerClass[: classNumber]) + componentNumber
                        variances[gaussianNumber, :, :] = variance
                        intervalSize = 2 * np.sqrt(np.diag(variance)) / numberOfComponents
                        means[gaussianNumber, :] = (mean - np.sqrt(np.diag(variance)) + intervalSize / 2 + (
                            componentNumber) * intervalSize).T
                        mixtureWeights[gaussianNumber] = 1 / numberOfComponents

            # Also remember the overall data variance for later usage in a conjugate prior on the variances
            dataMean = np.mean(data)
            tmp = data - dataMean
            dataVariance = np.var(tmp, axis=0)
            numberOfPseudoMeasurementsOfWishartPrior = 1
            pseudoVarianceOfWishartPrior = np.diag(dataVariance / numberOfPseudoMeasurementsOfWishartPrior)
            historyOfEMCost = [1 / eps]
            for EMIterationNumber in range(100):
                logger.debug('EMIterationNumber=%d', EMIterationNumber)

                # E-step: compute the posteriors based on the current parameters.

                for classNumber in range(numberOfClasses):
                    prior = priors[:, classNumber]
                    numberOfComponents = numberOfGaussiansPerClass[classNumber]
                    for componentNumber in range(numberOfComponents):
                        gaussianNumber = sum(numberOfGaussiansPerClass[:classNumber]) + componentNumber
                        mean = np.expand_dims(means[gaussianNumber, :], 1)
                        variance = variances[gaussianNumber, :, :]
                        L = np.linalg.cholesky(variance)
                        means_corrected_bias = biasCorrectedData.T - mean
                        if L.shape == (1, 1):
                            scale = 1.0 / L[0, 0]
                            tmp = means_corrected_bias * scale
                        else:
                            tmp = np.linalg.solve(L, means_corrected_bias)
                        tmp *= tmp
                        scaled_squared_mahalanobis_distances = np.sum(tmp, axis=0) * -0.5
                        sqrtDeterminantOfVariance = np.prod(np.diag(L))
                        scaling = 1.0 / (2 * np.pi) ** (
                                numberOfContrasts / 2) / sqrtDeterminantOfVariance
                        gaussianLikelihoods = np.exp(scaled_squared_mahalanobis_distances) * scaling
                        gaussianLikelihoods = gaussianLikelihoods.T
                        posteriors[:, gaussianNumber] = gaussianLikelihoods * (mixtureWeights[gaussianNumber] * prior)
                normalizer = np.sum(posteriors, axis=1) + eps
                posteriors = posteriors / np.expand_dims(normalizer, 1)

                minLogLikelihood = -np.sum(np.log(normalizer))
                intensityModelParameterCost = 0
                for gaussianNumber in range(numberOfGaussians):
                    variance = variances[gaussianNumber, :, :]
                    # Evaluate unnormalized Wishart distribution (conjugate prior on precisions) with parameters
                    #
                    #   scale matrix V = inv( pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior )
                    #
                    # and
                    #
                    #   degrees of freedom n = numberOfPseudoMeasurementsOfWishartPrior + numberOfContrasts + 1
                    #
                    # which has pseudoVarianceOfWishartPrior as the MAP solution in the absence of any data
                    #
                    minLogUnnormalizedWishart = np.trace(np.linalg.solve(variance,pseudoVarianceOfWishartPrior)) * \
                        numberOfPseudoMeasurementsOfWishartPrior / 2 + \
                        numberOfPseudoMeasurementsOfWishartPrior / 2 * np.log(np.linalg.det(variance))
                    intensityModelParameterCost = intensityModelParameterCost + minLogUnnormalizedWishart
                historyOfEMCost.append(minLogLikelihood + intensityModelParameterCost)

                priorEMCost = historyOfEMCost[-2]
                currentEMCost = historyOfEMCost[-1]
                costChangeEM = priorEMCost - currentEMCost
                changeCostEMPerVoxel = costChangeEM / activeVoxelCount
                changeCostEMPerVoxelThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion
                if changeCostEMPerVoxel < changeCostEMPerVoxelThreshold:
                    # Converged
                    print('EM converged!')
                    break

                # M-step: update the model parameters based on the current posterior
                #
                # First the mixture model parameters

                for gaussianNumber in range(numberOfGaussians):
                    posterior = posteriors[:, gaussianNumber]
                    posterior = posterior.reshape(-1, 1)
                    mean = biasCorrectedData.T @ posterior / np.sum(posterior)
                    tmp = biasCorrectedData - mean.T
                    variance = (tmp.T @ (tmp * posterior) + \
                                pseudoVarianceOfWishartPrior * numberOfPseudoMeasurementsOfWishartPrior) \
                               / (np.sum(posterior) + numberOfPseudoMeasurementsOfWishartPrior)
                    if modelSpecifications.useDiagonalCovarianceMatrices:
                        # Force diagonal covariance matrices
                        variance = np.diag(np.diag(variance))
                    variances[gaussianNumber, :, :] = variance
                    means[gaussianNumber, :] = mean.T
                mixtureWeights = np.sum(posteriors + eps, axis=0).T
                for classNumber in range(numberOfClasses):
                    # mixture weights are normalized (those belonging to one mixture sum to one)
                    numberOfComponents = numberOfGaussiansPerClass[classNumber]
                    gaussianNumbers = np.array(
                        np.sum(numberOfGaussiansPerClass[:classNumber]) + np.array(range(numberOfComponents)),
                        dtype=np.uint32)
                    mixtureWeights[gaussianNumbers] = mixtureWeights[gaussianNumbers] / np.sum(
                        mixtureWeights[gaussianNumbers])

                # Now update the parameters of the bias field model.
                if (estimateBiasField and (iterationNumber > 0)):

                    # Bias field correction: implements Eq. 8 in the paper
                    #    Van Leemput, "Automated Model-based Bias Field Correction of MR Images of the Brain", IEEE TMI 1999

                    precisions = np.zeros_like(variances)
                    for classNumber in range(numberOfGaussians):
                        precisions[classNumber, :, :] = np.linalg.inv(variances[classNumber, :, :]).reshape(
                            (1, numberOfContrasts, numberOfContrasts))
                    lhs = np.zeros((np.prod(numberOfBasisFunctions) * numberOfContrasts, np.prod(
                        numberOfBasisFunctions) * numberOfContrasts))  # left-hand side of linear system
                    rhs = np.zeros(
                        (np.prod(numberOfBasisFunctions) * numberOfContrasts, 1))  # right-hand side of linear system
                    weightsImageBuffer = np.zeros(downSampledImageSize)
                    tmpImageBuffer = np.zeros(downSampledImageSize)
                    numberOfBasisFunctions_prod = np.prod(numberOfBasisFunctions)
                    for contrastNumber1 in range(numberOfContrasts):
                        logger.debug('third time contrastNumber=%d', contrastNumber)
                        tmp = np.zeros((data.shape[0], 1), order='F')
                        for contrastNumber2 in range(numberOfContrasts):
                            classSpecificWeights = posteriors * precisions[:, contrastNumber1, contrastNumber2].T
                            weights = np.sum(classSpecificWeights, 1)
                            # Build up stuff needed for rhs
                            predicted = np.sum(classSpecificWeights * np.expand_dims(means[:, contrastNumber2], 2).T / (
                                    np.expand_dims(weights, 1) + eps), 1)
                            residue = data[:, contrastNumber2] - predicted
                            tmp = tmp + weights.reshape(-1, 1) * residue.reshape(-1, 1)
                            # Fill in submatrix of lhs
                            weightsImageBuffer[downSampledMaskIndices] = weights
                            computedPrecisionOfKroneckerProductBasisFunctions = computePrecisionOfKroneckerProductBasisFunctions(
                                downSampledKroneckerProductBasisFunctions, weightsImageBuffer)
                            lhs[
                            contrastNumber1 * numberOfBasisFunctions_prod: contrastNumber1 * numberOfBasisFunctions_prod + numberOfBasisFunctions_prod,
                            contrastNumber2 * numberOfBasisFunctions_prod:contrastNumber2 * numberOfBasisFunctions_prod + numberOfBasisFunctions_prod
                            ] = computedPrecisionOfKroneckerProductBasisFunctions
                        tmpImageBuffer[downSampledMaskIndices] = tmp.squeeze()
                        rhs[
                        contrastNumber1 * numberOfBasisFunctions_prod: contrastNumber1 * numberOfBasisFunctions_prod + numberOfBasisFunctions_prod] \
                            = projectKroneckerProductBasisFunctions(downSampledKroneckerProductBasisFunctions, tmpImageBuffer).reshape(-1, 1)
                    biasFieldCoefficients = np.linalg.solve(lhs, rhs).reshape((np.prod(numberOfBasisFunctions), numberOfContrasts), order='F')
                    downSampledBiasFields = biasCorrectData(biasCorrectedData, biasFieldCoefficients,
                                                            downSampledBiasCorrectedImageBuffers,
                                                            downSampledImageBuffers,
                                                            downSampledKroneckerProductBasisFunctions,
                                                            downSampledMask, downSampledMaskIndices,
                                                            numberOfContrasts)
                    pass
            if len(historyOfEMCost) > 2:
                visualizer.plot(historyOfEMCost[1:], title='History of EM Cost')
            visualizer.show( image_list=downSampledBiasFields, auto_scale=True, window_id='bias field', title='Bias Fields')
            visualizer.show( mesh=mesh, images=downSampledBiasCorrectedImageBuffers, window_id='samsegment em',
                title='Mesh Registration (EM)', names=[item.mergedName for item in modelSpecifications.sharedGMMParameters])
            historyOfEMCost = historyOfEMCost[1:]

            # Part II: update the position of the mesh nodes for the current mixture model and bias field parameter estimates

            downSampledBiasCorrectedImages = []
            for contrastNumber in range(numberOfContrasts):
                downSampledBiasCorrectedImages.append(gems.KvlImage(
                    requireNumpyArray(downSampledBiasCorrectedImageBuffers[:, :, :, contrastNumber])))

            # Set up cost calculator
            calculator = gems.KvlCostAndGradientCalculator(
                typeName='AtlasMeshToIntensityImage',
                images=downSampledBiasCorrectedImages,
                boundaryCondition='Sliding',
                transform=transform,
                means=means,
                variances=variances,
                mixtureWeights=mixtureWeights,
                numberOfGaussiansPerClass=numberOfGaussiansPerClass)

            optimizerType = 'L-BFGS'
            optimization_parameters = {
                'Verbose': optimizationOptions.verbose,
                'MaximalDeformationStopCriterion': optimizationOptions.maximalDeformationStopCriterion,
                'LineSearchMaximalDeformationIntervalStopCriterion': optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion,
                'MaximumNumberOfIterations': optimizationOptions.maximumNumberOfDeformationIterations,
                'BFGS-MaximumMemoryLength': optimizationOptions.BFGSMaximumMemoryLength
            }
            optimizer = gems.KvlOptimizer(optimizerType, mesh, calculator, optimization_parameters)
            historyOfDeformationCost = []
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            while True:
                minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
                print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                if maximalDeformation == 0:
                    break
                historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                historyOfMaximalDeformation.append(maximalDeformation)
            nodePositionsAfterDeformation = mesh.points
            maximalDeformationApplied = np.sqrt(
                np.max(np.sum((nodePositionsAfterDeformation - nodePositionsBeforeDeformation) ** 2, 1)))

            # print summary of iteration
            print('iterationNumber: %d' % iterationNumber)
            print('maximalDeformationApplied: %.4f' % maximalDeformationApplied)
            print('=======================================================')

            visualizer.show( mesh=mesh, images=downSampledBiasCorrectedImageBuffers, window_id='samsegment',
                title='Mesh Registration (Deformation)', names=[item.mergedName for item in modelSpecifications.sharedGMMParameters])

            # Keep track of the cost function we're optimizing
            historyOfCost.append(minLogLikelihoodTimesDeformationPrior + intensityModelParameterCost)
            priorCost = historyOfCost[-2]
            currentCost = historyOfCost[-1]
            costChange = priorCost - currentCost
            activeVoxelCount = len(downSampledMaskIndices[0])
            perVoxelDecrease = costChange / activeVoxelCount
            perVoxelDecreaseThreshold = optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion

            # Save history of the estimation
            if saveHistory:
                levelHistory['historyWithinEachIteration'].append({
                    'historyOfEMCost': historyOfEMCost,
                    'mixtureWeights': mixtureWeights,
                    'means': means,
                    'variances': variances,
                    'biasFieldCoefficients': biasFieldCoefficients,
                    'historyOfDeformationCost': historyOfDeformationCost,
                    'historyOfMaximalDeformation': historyOfMaximalDeformation,
                    'maximalDeformationApplied': maximalDeformationApplied
                })

            if perVoxelDecrease < perVoxelDecreaseThreshold:
                # Display the cost history
                if len(historyOfCost) > 2:
                    visualizer.plot(historyOfCost[1:], title='History of Cost')
                visualizer.show_movie(window_id='samsegment')
                # Log the final per-voxel cost
                with open(os.path.join(savePath, 'cost.txt'), "a") as file:
                    file.write("atlasRegistrationLevel%d %d %f\n" % (multiResolutionLevel, iterationNumber + 1, currentCost / activeVoxelCount))
                break

        # Get the final node positions
        finalNodePositions = mesh.points

        # Transform back in template space (i.e., undoing the affine registration that we applied), and save for later usage
        tmp = np.linalg.solve(totalTransformationMatrix, np.pad(finalNodePositions, ((0, 0), (0, 1)), 'constant', constant_values=1).T).T
        finalNodePositionsInTemplateSpace = tmp[:, 0: 3]

        # Record deformation delta here in lieu of maintaining history
        nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel = finalNodePositionsInTemplateSpace - initialNodePositionsInTemplateSpace

        # Save history of the estimation
        if saveHistory:
            levelHistory['downSamplingFactors'] = downSamplingFactors
            levelHistory['downSampledImageBuffers'] = downSampledImageBuffers
            levelHistory['downSampledMask'] = downSampledMask
            levelHistory['initialNodePositions'] = initialNodePositions
            levelHistory['finalNodePositions'] = finalNodePositions
            levelHistory['initialNodePositionsInTemplateSpace'] = initialNodePositionsInTemplateSpace
            levelHistory['finalNodePositionsInTemplateSpace'] = finalNodePositionsInTemplateSpace
            levelHistory['historyOfCost'] = historyOfCost
            levelHistory['priorsAtEnd'] = priors
            levelHistory['posteriorsAtEnd'] = posteriors
            history['historyWithinEachMultiResolutionLevel'].append(levelHistory)

    # End resolution level loop

    # Save the history
    if saveHistory:
        with open(os.path.join(savePath, 'history.p'), 'wb') as file:
            pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)

    # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
    # with all the original labels (instead of the reduced "super"-structure labels we created)

    fs.printPeakMemory('samsegment starting segmentation')

    # Get bias field corrected images
    biasCorrectedImageBuffers = np.zeros((imageSize[0], imageSize[1], imageSize[2], numberOfContrasts))
    biasFields = np.zeros((imageSize[0], imageSize[1], imageSize[2], numberOfContrasts))

    for contrastNumber in range(numberOfContrasts):
        biasField = backprojectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, biasFieldCoefficients[:, contrastNumber])
        biasCorrectedImageBuffers[:, :, :, contrastNumber] = imageBuffers[:, :, :, contrastNumber] - biasField * mask
        biasFields[:, :, :, contrastNumber] = biasField

    # Read the atlas, applying the affine registration transform
    mesh_collection = gems.KvlMeshCollection()
    mesh_collection.read(modelSpecifications.atlasFileName)
    mesh_collection.k = modelSpecifications.K
    mesh_collection.transform(transform)
    mesh = mesh_collection.reference_mesh

    # Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
    nodePositions = mesh.points
    numberOfNodes = nodePositions.shape[0]
    transformMatrix = transform.as_numpy_array
    tmp = np.linalg.solve(transformMatrix, np.pad(nodePositions, ((0, 0), (0, 1)), mode='constant', constant_values=1).T).T
    nodePositionsInTemplateSpace = tmp[:, 0: 3]

    # Get the estimated warp in template space
    [estimatedNodeDeformationInTemplateSpace, estimated_averageDistance, estimated_maximumDistance] = gems.kvlWarpMesh(
        optimizationOptions.multiResolutionSpecification[-1].atlasFileName,
        nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel,
        modelSpecifications.atlasFileName
    )

    # Apply this warp on the mesh node positions in template space, and transform into current space
    desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace
    tmp = (transformMatrix @ np.pad(desiredNodePositionsInTemplateSpace, ((0, 0), (0, 1)), mode='constant', constant_values=1).T).T
    desiredNodePositions = tmp[:, 0: 3]
    mesh.points = requireNumpyArray(desiredNodePositions)
    alphas = mesh.alphas
    numberOfStructures = alphas.shape[1]

    # Get the priors as dictated by the current mesh position
    data = biasCorrectedImageBuffers
    priors = mesh.rasterize_3(imageSize, -1)

    # NOTE: NOT GOING TO RESHAPE, WILL USE MASK INDEXING
    # Ignore everything that's has zero intensity
    priors = priors[mask, :]
    data = data[mask, :]
    likelihood_count = data.shape[0]

    # Calculate the posteriors
    posteriors = np.zeros_like(priors, dtype=np.float64)
    for structureNumber in range(numberOfStructures):
        prior = priors[:, structureNumber] / 65535
        mixedLikelihoods = np.zeros((likelihood_count, 1))
        for classNumber in range(numberOfClasses):
            fraction = translationTable[classNumber, structureNumber]
            if fraction < 1e-10: continue
            # Compute likelihood of this class (aka mixture model)
            likelihoods = np.zeros((likelihood_count, 1))
            numberOfComponents = numberOfGaussiansPerClass[classNumber]
            for componentNumber in range(numberOfComponents):
                gaussianNumber = int(np.sum(numberOfGaussiansPerClass[: classNumber]) + componentNumber)
                mean = np.expand_dims(ensureDims(means, 2)[gaussianNumber, :], 1)
                variance = ensureDims(variances, 3)[gaussianNumber, :, :]
                mixtureWeight = mixtureWeights[gaussianNumber]
                L = np.linalg.cholesky(variance)
                tmp = np.linalg.solve(L, data.T - mean)
                squaredMahalanobisDistances = (np.sum(tmp ** 2, axis=0)).T
                sqrtDeterminantOfVariance = np.prod(np.diag(L))
                gaussianLikelihoods = np.exp(-squaredMahalanobisDistances / 2) / (2 * np.pi) ** (
                        numberOfContrasts / 2) / sqrtDeterminantOfVariance
                likelihoods = likelihoods + ensureDims(gaussianLikelihoods, 2) * mixtureWeight
            mixedLikelihoods = mixedLikelihoods + likelihoods * fraction
        posteriors[:, structureNumber] = np.squeeze(mixedLikelihoods) * prior
    normalizer = np.sum(posteriors, 1) + eps
    posteriors = posteriors / ensureDims(normalizer, 2)

    # Compute volumes in mm^3
    volumeOfOneVoxel = np.abs(np.linalg.det(imageToWorldTransformMatrix[0:3, 0:3]))
    volumesInCubicMm = (np.sum(posteriors, axis=0)) * volumeOfOneVoxel

    # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
    structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)
    freeSurferSegmentation = np.zeros(imageSize, dtype=np.uint16)
    FreeSurferLabels = np.array(FreeSurferLabels, dtype=np.uint16)
    freeSurferSegmentation[mask] = FreeSurferLabels[structureNumbers]

    # Write to file, remembering to un-crop the segmentation to the original image size
    uncroppedFreeSurferSegmentation = np.zeros(nonCroppedImageSize, dtype=np.float32)
    uncroppedFreeSurferSegmentation[cropping] = freeSurferSegmentation
    print('Writing out freesurfer segmentation')
    gems.KvlImage(requireNumpyArray(uncroppedFreeSurferSegmentation)).write(
        os.path.join(savePath, 'crispSegmentation.nii'),
        gems.KvlTransform(requireNumpyArray(imageToWorldTransformMatrix))
    )

    # Save the final mesh collection
    if saveMesh:
        print('Saving the final mesh in template space')
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(modelSpecifications.atlasFileName)
        mesh_collection.reference_mesh.points = requireNumpyArray(desiredNodePositionsInTemplateSpace)
        mesh_collection.write(os.path.join(savePath, 'mesh.txt'))

    # Isolate the white matter posteriors to compute the average intensity of white matter voxels for each input
    wmRightIndex = names.index('Right-Cerebral-White-Matter')
    wmLeftIndex = names.index('Left-Cerebral-White-Matter')
    wmPosteriors = np.zeros(imageSize, dtype=float)
    wmPosteriors[mask] = posteriors[:, wmRightIndex] + posteriors[:, wmLeftIndex]

    # Also write out the bias field and the bias corrected image, each time remembering to un-crop the images
    for contrastNumber, imageFileName in enumerate(imageFileNames):
        image_base_path, ext = os.path.splitext(imageFileName)
        data_path, scanName = os.path.split(image_base_path)

        # First bias field - we're computing it also outside of the mask, but clip the intensities there to
        # the range observed inside the mask (with some margin) to avoid crazy extrapolation values
        logBiasField = biasFields[:, :, :, contrastNumber]
        clippingMargin = np.log(2)
        clippingMin = logBiasField[mask].min() - clippingMargin
        clippingMax = logBiasField[mask].max() + clippingMargin
        logBiasField[ logBiasField < clippingMin ] = clippingMin
        logBiasField[ logBiasField > clippingMax ] = clippingMax
        biasField = np.zeros(nonCroppedImageSize, dtype=np.float32)
        biasField[cropping] = np.exp( logBiasField )
        outputFileName = os.path.join(savePath, scanName + '_biasField.nii')
        gems.KvlImage(biasField).write(outputFileName, gems.KvlTransform(requireNumpyArray(imageToWorldTransformMatrix)))

        # Then save the bias-corrected image (scale the white matter average to 110)
        biasCorrected = np.exp(biasCorrectedImageBuffers[:, :, :, contrastNumber])
        biasCorrectedOutput = np.zeros(nonCroppedImageSize, dtype=np.float32)
        scaleFactor = 110 / np.average(biasCorrected, weights=wmPosteriors)
        biasCorrectedOutput[cropping] = biasCorrected * scaleFactor
        outputFileName = os.path.join(savePath, scanName + '_biasCorrected.nii')
        gems.KvlImage(biasCorrectedOutput).write(outputFileName, gems.KvlTransform(requireNumpyArray(imageToWorldTransformMatrix)))

        # Save a note indicating the scaling factor
        with open(os.path.join(savePath, 'scaling-factor.txt'), 'w') as f:
            print(scaleFactor, file=f)

    return [FreeSurferLabels, names, volumesInCubicMm]
