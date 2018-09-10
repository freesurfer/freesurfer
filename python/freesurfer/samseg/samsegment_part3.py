import os
import numpy as np
import freesurfer.gems as gems

from .utilities import require_np_array
from .bias_correction import backprojectKroneckerProductBasisFunctions
from .figures import initVisualizer

# from .dev_utils.debug_client import run_test_cases, create_checkpoint_manager, create_part3_inspection_team, load_starting_fixture


eps = np.finfo(float).eps


def ensure_dims(np_array, dims):
    if np_array.ndim < dims:
        return ensure_dims(np.expand_dims(np_array, axis=dims), dims)
    elif np_array.ndim == dims:
        return np_array


def samsegment_part3(
        modelSpecifications,
        optimizationOptions,
        part1_results_dict,
        part2_results_dict,
        imageFileNames,
        visualizer,
        checkpoint_manager=None
):

    # Setup null visualization if necessary
    if visualizer is None: visualizer = initVisualizer(False, False)

    croppingOffset = part1_results_dict['croppingOffset']
    FreeSurferLabels = part1_results_dict['FreeSurferLabels']
    imageSize = part1_results_dict['imageSize']
    imageToWorldTransformMatrix = part1_results_dict['imageToWorldTransformMatrix']
    kroneckerProductBasisFunctions = part1_results_dict['kroneckerProductBasisFunctions']
    mask = part1_results_dict['mask']
    nonCroppedImageSize = part1_results_dict['nonCroppedImageSize']
    numberOfClasses = part1_results_dict['numberOfClasses']
    numberOfContrasts = part1_results_dict['numberOfContrasts']
    numberOfGaussiansPerClass = part1_results_dict['numberOfGaussiansPerClass']
    translationTable = part1_results_dict['translationTable']
    savePath = part1_results_dict['savePath']

    biasFieldCoefficients = part2_results_dict['biasFieldCoefficients']
    imageBuffers = part2_results_dict['imageBuffers']
    means = part2_results_dict['means']
    mixtureWeights = part2_results_dict['mixtureWeights']
    nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel = \
        part2_results_dict['nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel']
    variances = part2_results_dict['variances']
    transform = gems.KvlTransform(np.asfortranarray(part2_results_dict['transformMatrix']))

    nonCroppedImageSize = [int(dim) for dim in nonCroppedImageSize]
    croppingOffset = [int(offset) for offset in croppingOffset]

    # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
    # with all the original labels (instead of the reduced "super"-structure labels we created).
    #
    # Get bias field corrected images
    biasCorrectedImageBuffers = np.zeros((imageSize[0], imageSize[1], imageSize[2], numberOfContrasts))
    biasFields = np.zeros((imageSize[0], imageSize[1], imageSize[2], numberOfContrasts))

    for contrastNumber in range(numberOfContrasts):
        biasField = backprojectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions,
                                                              biasFieldCoefficients[:, contrastNumber])
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
    tmp = np.linalg.solve(transformMatrix,
                          np.pad(nodePositions, ((0, 0), (0, 1)), mode='constant', constant_values=1).T).T
    nodePositionsInTemplateSpace = tmp[:, 0: 3]
    
    # Get the estimated warp in template space
    [estimatedNodeDeformationInTemplateSpace, estimated_averageDistance, estimated_maximumDistance] = gems.kvlWarpMesh(
        optimizationOptions.multiResolutionSpecification[-1].atlasFileName,
        nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel,
        modelSpecifications.atlasFileName
    )

    # Apply this warp on the mesh node positions in template space, and transform into current space
    desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace
    tmp = (transformMatrix @ np.pad(desiredNodePositionsInTemplateSpace, ((0, 0), (0, 1)), mode='constant',
                                    constant_values=1).T).T
    desiredNodePositions = tmp[:, 0: 3]
    mesh.points = require_np_array(desiredNodePositions)
    alphas = mesh.alphas
    numberOfStructures = alphas.shape[1]
    
    # Get the priors as dictated by the current mesh position
    data = biasCorrectedImageBuffers
    priors = mesh.rasterize_3(imageSize, -1)
    
    # NOT GOING TO RESHAPE, WILL USE MASK INDEXING
    #
    #
    # Ignore everything that's has zero intensity
    priors = priors[mask == 1, :]
    data = data[mask == 1, :]
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
                mean = np.expand_dims(ensure_dims(means, 2)[gaussianNumber, :], 1)
                variance = ensure_dims(variances, 3)[gaussianNumber, :, :]
                mixtureWeight = mixtureWeights[gaussianNumber]
                L = np.linalg.cholesky(variance)
                tmp = np.linalg.solve(L, data.T - mean)
                squaredMahalanobisDistances = (np.sum(tmp ** 2, axis=0)).T
                sqrtDeterminantOfVariance = np.prod(np.diag(L))
                gaussianLikelihoods = np.exp(-squaredMahalanobisDistances / 2) / (2 * np.pi) ** (
                        numberOfContrasts / 2) / sqrtDeterminantOfVariance
                likelihoods = likelihoods + ensure_dims(gaussianLikelihoods, 2) * mixtureWeight
            mixedLikelihoods = mixedLikelihoods + likelihoods * fraction
        posteriors[:, structureNumber] = np.squeeze(mixedLikelihoods) * prior
    normalizer = np.sum(posteriors, 1) + eps
    posteriors = posteriors / ensure_dims(normalizer, 2)
    
    # Compute volumes in mm^3
    volumeOfOneVoxel = np.abs(np.linalg.det(imageToWorldTransformMatrix[0:3, 0:3]))
    volumesInCubicMm = (np.sum(posteriors, axis=0)) * volumeOfOneVoxel
    
    # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
    structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)
    freeSurferSegmentation = np.zeros(imageSize, dtype=np.uint16)
    FreeSurferLabels = np.array(FreeSurferLabels, dtype=np.uint16)
    freeSurferSegmentation[mask == 1] = FreeSurferLabels[structureNumbers]

    # Write to file, remembering to un-crop the segmentation to the original image size
    uncroppedFreeSurferSegmentation = np.zeros(nonCroppedImageSize, dtype=np.float32)
    uncroppedFreeSurferSegmentation[croppingOffset[0]: imageSize[0] + croppingOffset[0],
    croppingOffset[1]: imageSize[1] + croppingOffset[1],
    croppingOffset[2]: imageSize[2] + croppingOffset[2]] = freeSurferSegmentation
    print('Writing out freesurfer segmentation')
    gems.KvlImage(require_np_array(uncroppedFreeSurferSegmentation)).write(
        os.path.join(savePath, 'crispSegmentation.nii'),
        gems.KvlTransform(require_np_array(imageToWorldTransformMatrix))
    )
    
    # Also write out the bias field and the bias corrected image, each time remembering to un-crop the images
    for contrastNumber, imageFileName in enumerate(imageFileNames):
        image_base_path, ext = os.path.splitext(imageFileName)
        data_path, scanName = os.path.split(image_base_path)
        
        #  First bias field
        biasField = np.zeros(nonCroppedImageSize, dtype=np.float32)
        biasField[
            croppingOffset[0]:croppingOffset[0] + imageSize[0],
            croppingOffset[1]:croppingOffset[1] + imageSize[1],
            croppingOffset[2]:croppingOffset[2] + imageSize[2],
        ] = np.exp(biasFields[:, :, :, contrastNumber]) * mask
        outputFileName = os.path.join(savePath, scanName + '_biasField.nii')
        gems.KvlImage(biasField).write(
            outputFileName,
            gems.KvlTransform(require_np_array(imageToWorldTransformMatrix))
        )
        
        #  Then bias field corrected data
        biasCorrected = np.zeros(nonCroppedImageSize, dtype=np.float32)
        biasCorrected[
            croppingOffset[0]:croppingOffset[0] + imageSize[0],
            croppingOffset[1]:croppingOffset[1] + imageSize[1],
            croppingOffset[2]:croppingOffset[2] + imageSize[2],
        ] = np.exp(biasCorrectedImageBuffers[:, :, :, contrastNumber])
        outputFileName = os.path.join(savePath, scanName + '_biasCorrected.nii')
        gems.KvlImage(biasCorrected).write(
            outputFileName,
            gems.KvlTransform(require_np_array(imageToWorldTransformMatrix))
        )
    return {
        'FreeSurferLabels': FreeSurferLabels,
        'freeSurferSegmentation': freeSurferSegmentation,
        'uncroppedFreeSurferSegmentation': uncroppedFreeSurferSegmentation,
        'volumesInCubicMm': volumesInCubicMm,
    }


def test_samseg_ported_part3(case_name, case_file_folder, savePath):
    checkpoint_manager = create_checkpoint_manager(case_file_folder)
    fixture = load_starting_fixture()
    part1_results_dict, part1_results_dict_python, part1_results_dict_matlab = checkpoint_manager.substitute('part1', 1)
    part1_results_dict['savePath'] = savePath
    part2_results_dict, part2_results_dict_python, part2_results_dict_matlab = checkpoint_manager.substitute('part2', 1)
    part3_results_dict = samsegment_part3(
        fixture['modelSpecifications'],
        fixture['optimizationOptions'],
        part1_results_dict,
        part2_results_dict,
        fixture['imageFileNames'],
        checkpoint_manager
    )
    checkpoint_manager.save(part3_results_dict, 'part3', 1)
    names = part1_results_dict['names']
    FreeSurferLabels = part3_results_dict['FreeSurferLabels']
    volumesInCubicMm = part3_results_dict['volumesInCubicMm']
    print('names', names)
    print('FreeSurferLabels', FreeSurferLabels)
    print('volumesInCubicMm', volumesInCubicMm)
    create_part3_inspection_team().inspect_all(checkpoint_manager)


if __name__ == '__main__':
    run_test_cases(action=test_samseg_ported_part3)
