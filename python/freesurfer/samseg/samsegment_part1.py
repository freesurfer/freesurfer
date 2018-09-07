import math
import numpy as np
from functools import reduce
from operator import mul

import freesurfer.gems as gems

from .figures import initVisualizer

# from .dev_utils.debug_client import create_checkpoint_manager, run_test_cases, create_part1_inspection_team, load_starting_fixture


def samsegment_part1(
        imageFileNames,
        transformedTemplateFileName,
        modelSpecifications,
        optimizationOptions,
        savePath,
        visualizer,
        checkpoint_manager=None
):
    
    # Setup null visualization if necessary
    if visualizer is None: visualizer = initVisualizer(False, False)

    # Print input options
    print('==========================')
    print_image_file_names(imageFileNames)
    print_transformed_template_file_name(transformedTemplateFileName)
    print_model_specifications(modelSpecifications)
    print_optimization_options(optimizationOptions)
    print('savePath:', savePath)

    numberOfContrasts = len(imageFileNames)
    # Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
    # translation, rotation, scaling, and skewing) as well - this transformation will later be used
    # to initially transform the location of the atlas mesh's nodes into the coordinate system of
    # the image.

    imageBuffers = []
    images = []
    for imageFileName in imageFileNames:
        # Get the pointers to image and the corresponding transform
        image = gems.KvlImage(imageFileName, transformedTemplateFileName)
        transform = image.transform_matrix
        nonCroppedImageSize = image.non_cropped_image_size
        croppingOffset = image.cropping_offset
        images.append(image)
        imageBuffers.append(image.getImageBuffer())
    nonCroppedImageSize = [int(dim) for dim in nonCroppedImageSize]
    croppingOffset = [int(offset) for offset in croppingOffset]

    imageSize = imageBuffers[0].shape
    imageBuffers = np.transpose(imageBuffers, axes=[1, 2, 3, 0])
    visualizer.show(images=imageBuffers, window_id='samsegment contrast', title='Samsegment Contrasts')

    # Also read in the voxel spacing -- this is needed since we'll be specifying bias field smoothing kernels, downsampling
    # steps etc in mm.
    nonCroppedImage = gems.KvlImage(imageFileNames[0])
    imageToWorldTransformMatrix = nonCroppedImage.transform_matrix.as_numpy_array
    voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)

    # Read the atlas mesh from file, immediately applying the previously determined transform to the location
    # of its nodes. Rather than one mesh, the atlas consists of a so-called "collection" of meshes: they're
    # all the same except for the location of their mesh nodes. The first mesh, so-called "reference mesh"
    # has index -1: it represents the average shape (location of the mesh nodes) of all the other meshes. If
    # the atlas was built from N training subjects, there will be N other meshes as well, each one warping
    # the atlas optimally to each of the training subjects
    #
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

    visualizer.show(
        probabilities=backgroundPrior,
        images=imageBuffers,
        window_id='samsegment background',
        title='Samsegment Background Priors'
    )
    smoothingSigmas = [1.0 * modelSpecifications.brainMaskingSmoothingSigma] * 3
    smoothedBackgroundPrior = gems.KvlImage.smooth_image_buffer(backgroundPrior, smoothingSigmas)
    visualizer.show(
        probabilities=smoothedBackgroundPrior,
        window_id='samsegment smoothed',
        title='Samsegment Smoothed Background Priors'
    )

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
    for contrastNumber in range(numberOfContrasts):
        imageBuffers[:, :, :, contrastNumber] *= brainMask

    visualizer.show(
        images=imageBuffers,
        window_id='samsegment images',
        title='Samsegment Masked Contrasts'
    )

    # Let's prepare for the bias field correction that is part of the imaging model. It assumes
    # an additive effect, whereas the MR physics indicate it's a multiplicative one - so we log
    # transform the data first. In order to do so, mask out zeros from
    # the images.
    # This removes any voxel where any contrast has a zero value
    # (messes up log)
    mask = np.full(imageSize, True, dtype=np.bool)
    for contrastNumber in range(numberOfContrasts):
        ## Note maskIndices not needed for python port.
        mask = mask * (imageBuffers[:, :, :, contrastNumber] > 0)
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings('ignore')
        log_buffers = np.log(imageBuffers)

    imageBuffers = np.ma.fix_invalid(log_buffers).filled(0)
    for contrastNumber in range(numberOfContrasts):
        imageBuffers[np.logical_not(mask), contrastNumber] = 0
    log_buffers = None

    # Merge classes into "super-structures" that define a single Gaussian mixture model shared between the classes belonging
    # to the same super-structure
    FreeSurferLabels = modelSpecifications.FreeSurferLabels
    names = modelSpecifications.names
    colors = modelSpecifications.colors
    [reducedAlphas, reducedNames, reducedFreeSurferLabels, reducedColors, translationTable
     ] = gems.kvlMergeAlphas(alphas, names, modelSpecifications.sharedGMMParameters, FreeSurferLabels, colors)

    visualizer.show(
        mesh=mesh,
        shape=imageBuffers.shape,
        window_id='samsegment mesh',
        title='Samsegment Mesh',
        names=names,
        legend_width=350, # Wider legend to handle the longer names
    )

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
    #
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
    #
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
    return {
        'biasFieldCoefficients': biasFieldCoefficients,
        'colors': colors,
        'croppingOffset': croppingOffset,
        'FreeSurferLabels': FreeSurferLabels,
        'imageBuffers': imageBuffers,
        'imageSize': imageSize,
        'imageToWorldTransformMatrix': imageToWorldTransformMatrix,
        'kroneckerProductBasisFunctions': kroneckerProductBasisFunctions,
        'mask': mask,
        'names': names,
        'nonCroppedImageSize': nonCroppedImageSize,
        'numberOfBasisFunctions': numberOfBasisFunctions,
        'numberOfClasses': numberOfClasses,
        'numberOfContrasts': numberOfContrasts,
        'numberOfGaussians': numberOfGaussians,
        'numberOfGaussiansPerClass': numberOfGaussiansPerClass,
        'translationTable': translationTable,
        'savePath': savePath,
        'transformMatrix': transform.as_numpy_array,
        'voxelSpacing': voxelSpacing,
    }


def print_image_file_names(imageFileNames):
    print('imageFileNames: ')
    for imageFileName in imageFileNames:
        print('    '.format(imageFileName))
    print('-----')


def print_transformed_template_file_name(transformedTemplateFileName):
    print('transformedTemplateFileName: {0}'.format(transformedTemplateFileName))
    print()
    print('-----')


def print_model_specifications(modelSpecifications):
    print('modelSpecifications: ')
    print('    atlasFileName={0}'.format(modelSpecifications.atlasFileName))
    print('    useDiagonalCovarianceMatrices={0}'.format(modelSpecifications.useDiagonalCovarianceMatrices))
    print('    brainMaskingSmoothingSigma={0}'.format(modelSpecifications.brainMaskingSmoothingSigma))
    print('    brainMaskingThreshold={0}'.format(modelSpecifications.brainMaskingThreshold))
    print('    K={0}'.format(modelSpecifications.K))
    print('    biasFieldSmoothingKernelSize={0}'.format(modelSpecifications.biasFieldSmoothingKernelSize))
    print('    sharedGMMParameters:')
    for params in modelSpecifications.sharedGMMParameters:
        print('        mergedName={0}'.format(params.mergedName))
        print('        number of components={0}'.format(params.numberOfComponents))
        print('        searchStrings={0}'.format(params.searchStrings))
        print('        ---')
    print('-----')


def print_optimization_options(options):
    print('optimizationOptions:')
    print('    maximumNumberOfDeformationIterations={0}'.format(options.maximumNumberOfDeformationIterations))
    print('    absoluteCostPerVoxelDecreaseStopCriterion={0}'.format(options.absoluteCostPerVoxelDecreaseStopCriterion))
    print('    verbose={0}'.format(options.verbose))
    print('    maximalDeformationStopCriterion={0}'.format(options.maximalDeformationStopCriterion))
    print('    lineSearchMaximalDeformationIntervalStopCriterion={0}'.format(
        options.lineSearchMaximalDeformationIntervalStopCriterion))
    print('    maximalDeformationAppliedStopCriterion={0}'.format(options.maximalDeformationAppliedStopCriterion))
    print('    BFGSMaximumMemoryLength={0}'.format(options.BFGSMaximumMemoryLength))
    print('    multiResolutionSpecification:')
    for spec in options.multiResolutionSpecification:
        print('        atlasFileName={0}'.format(spec.atlasFileName))
        print('        targetDownsampledVoxelSpacing={0}'.format(spec.targetDownsampledVoxelSpacing))
        print('        maximumNumberOfIterations={0}'.format(spec.maximumNumberOfIterations))
        print('        estimateBiasField={0}'.format(spec.estimateBiasField))
        print('        ---')
    print('-----')


def test_samseg_part_1(case_name, case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    fixture = load_starting_fixture()
    part0_results_dict = checkpoint_manager.load_specification('part0', 1)
    part0_results_dict, part0_results_dict_python, part0_results_dict_matlab = checkpoint_manager.substitute('part0', 1,
                                                                                                             python_dict=part0_results_dict)
    part0_results_dict['savePath'] = savePath

    part1_results_dict = samsegment_part1(
        part0_results_dict['imageFileNames'],
        part0_results_dict['transformedTemplateFileName'],
        fixture['modelSpecifications'],
        fixture['optimizationOptions'],
        savePath,
        fixture['visualizer'],
        checkpoint_manager
    )
    create_part1_inspection_team().inspect_all(checkpoint_manager)


if __name__ == '__main__':
    run_test_cases(action=test_samseg_part_1)
