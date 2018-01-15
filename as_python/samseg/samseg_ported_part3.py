import os

import GEMS2Python
import numpy as np

from as_python.samseg.bias_correction import backprojectKroneckerProductBasisFunctions
from as_python.samseg.dev_utils.debug_client import CheckpointManager
from as_python.samseg.kvlWarpMesh import kvlWarpMesh

# from dotdict import DotDict
from as_python.samseg.run_utilities import load_starting_fixture

eps = np.finfo(float).eps


def ensure_dims(np_array, dims):
    if np_array.ndim < dims:
        return ensure_dims(np.expand_dims(np_array, axis=dims), dims)
    elif np_array.ndim == dims:
        return np_array

def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])

def samsegment_part3(
    modelSpecifications,
    optimizationOptions,
    part1_results_dict,
    part2_results_dict,
    checkpoint_manager=None
):
    biasFieldCoefficients = part1_results_dict['biasFieldCoefficients']
    croppingOffset = part1_results_dict['croppingOffset']
    FreeSurferLabels = part1_results_dict['FreeSurferLabels']
    imageSize = part1_results_dict['imageSize']
    imageToWorldTransformMatrix = part1_results_dict['imageToWorldTransformMatrix']
    kroneckerProductBasisFunctions = part1_results_dict['kroneckerProductBasisFunctions']
    mask = part1_results_dict['mask']
    nonCroppedImageSize = part1_results_dict['nonCroppedImageSize']
    numberOfContrasts = part1_results_dict['numberOfContrasts']
    numberOfGaussiansPerClass = part1_results_dict['numberOfGaussiansPerClass']
    reducingLookupTable = part1_results_dict['reducingLookupTable']
    savePath = part1_results_dict['savePath']

    imageBuffers = part2_results_dict['imageBuffers']
    means = part2_results_dict['means']
    mixtureWeights = part2_results_dict['mixtureWeights']
    nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel = \
        part2_results_dict['nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel']
    variances = part2_results_dict['variances']
    transform = GEMS2Python.KvlTransform(np.asfortranarray(part2_results_dict['transformMatrix']))

    nonCroppedImageSize = [int(dim) for dim in nonCroppedImageSize]
    croppingOffset = [int(offset) for offset in croppingOffset]

    # ﻿
    # % Save something about how the estimation proceeded
    # history.imageBuffers = imageBuffers;
    # history.mask = mask;
    # history.historyWithinEachMultiResolutionLevel = historyWithinEachMultiResolutionLevel;
    # eval( [ 'save ' savePath '/history.mat history -v7.3' ] );
    #
    #
    # % OK, now that all the parameters have been estimated, try to segment the original, full resolution image
    # % with all the original labels (instead of the reduced "super"-structure labels we created).
    #
    # % Get bias field corrected images
    # biasCorrectedImageBuffers = zeros( [ imageSize numberOfContrasts ] );

    biasCorrectedImageBuffers = np.zeros( (imageSize[0], imageSize[1], imageSize[2], numberOfContrasts) )
    # biasFields = zeros( [ imageSize numberOfContrasts ] );
    biasFields = np.zeros( (imageSize[0], imageSize[1], imageSize[2], numberOfContrasts) )

    # TODO remove these ensure_dims once merging with part 2
    biasFieldCoefficients = ensure_dims(biasFieldCoefficients, 2)
    imageBuffers = ensure_dims(imageBuffers, 4)
    # for contrastNumber = 1 : numberOfContrasts
    for contrastNumber in range(numberOfContrasts):
        #   biasField = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
        biasField = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, biasFieldCoefficients[ :, contrastNumber ] )
        #   biasCorrectedImageBuffers( :, :, :, contrastNumber ) = imageBuffers( :, :, :, contrastNumber ) - biasField .* mask;
        biasCorrectedImageBuffers[ :, :, :, contrastNumber ] = imageBuffers[:, :, :, contrastNumber] - biasField * mask
        #   biasFields( :, :, :, contrastNumber ) = biasField;
        biasFields[:, :, :, contrastNumber] = biasField
        # end
    #
    #
    #
    # % Read the atlas, applying the affine registration transform
    # meshCollection = kvlReadMeshCollection( modelSpecifications.atlasFileName, transform, modelSpecifications.K );

    mesh_collection = GEMS2Python.KvlMeshCollection()
    mesh_collection.read(modelSpecifications.atlasFileName)
    mesh_collection.k = modelSpecifications.K
    mesh_collection.transform(transform)


    #   mesh = kvlGetMesh( meshCollection, -1 );
    mesh = mesh_collection.reference_mesh
    #
    # % Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
    # nodePositions = kvlGetMeshNodePositions( mesh );
    nodePositions = mesh.points
    # numberOfNodes = size( nodePositions, 1 );
    numberOfNodes = nodePositions.shape[0]
    # transformMatrix = double( kvlGetTransformMatrix( transform ) );
    transformMatrix = transform.as_numpy_array
    # tmp = ( transformMatrix \ [ nodePositions ones( numberOfNodes, 1 ) ]' )';
    tmp = np.linalg.solve(transformMatrix, np.pad(nodePositions, ((0,0), (0,1)),  mode='constant', constant_values=1).T).T
    # nodePositionsInTemplateSpace = tmp( :, 1 : 3 );
    nodePositionsInTemplateSpace = tmp[:, 0 : 3 ]
    #
    # % Get the estimated warp in template space
    # estimatedNodeDeformationInTemplateSpace = ...
    #           kvlWarpMesh( optimizationOptions.multiResolutionSpecification( end ).atlasFileName, ...
    #                         historyWithinEachMultiResolutionLevel( end ).finalNodePositionsInTemplateSpace ...
    #                         - historyWithinEachMultiResolutionLevel( end ).initialNodePositionsInTemplateSpace, ...
    #                         modelSpecifications.atlasFileName );
    #
    [estimatedNodeDeformationInTemplateSpace, estimated_averageDistance, estimated_maximumDistance] = kvlWarpMesh(
        optimizationOptions.multiResolutionSpecification[-1].atlasFileName,
        nodeDeformationInTemplateSpaceAtPreviousMultiResolutionLevel,
        modelSpecifications.atlasFileName
    )

    # % Apply this warp on the mesh node positions in template space, and transform into current space
    # desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace;
    desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace
    # tmp = ( transformMatrix * [ desiredNodePositionsInTemplateSpace ones( numberOfNodes, 1 ) ]' )';
    tmp = ( transformMatrix @ np.pad(desiredNodePositionsInTemplateSpace, ((0,0), (0,1)),  mode='constant', constant_values=1).T).T
    # desiredNodePositions = tmp( :, 1 : 3 );
    desiredNodePositions = tmp[ :, 0 : 3 ]
    #
    # %
    # kvlSetMeshNodePositions( mesh, desiredNodePositions );
    mesh.points = require_np_array(desiredNodePositions)
    #
    # %
    # alphas = kvlGetAlphasInMeshNodes( mesh );
    alphas = mesh.alphas
    # numberOfStructures = size( alphas, 2 );
    numberOfStructures = alphas.shape[1]
    #
    #
    #
    # % Get the priors as dictated by the current mesh position
    # data = reshape( biasCorrectedImageBuffers, [ prod( imageSize ) numberOfContrasts ] );
    data = biasCorrectedImageBuffers
    # priors = kvlRasterizeAtlasMesh( mesh, imageSize );
    priors = mesh.rasterize(imageSize, -1)
    # priors = reshape( priors, [ prod( imageSize ) numberOfStructures ] );
    # NOT GOING TO RESHAPE, WILL USE MASK INDEXING
    #
    #
    # % Ignore everything that's has zero intensity
    # priors = priors( maskIndices, : );
    priors = priors[mask == 1, :]
    # data = data( maskIndices, : );
    data = data[mask == 1, :]
    likelihood_count = data.shape[0]
    #
    #
    # % Calculate the posteriors
    # posteriors = zeros( size( priors ), 'double' );
    posteriors = np.zeros_like(priors, dtype=np.float64)
    # for structureNumber = 1 : numberOfStructures
    for structureNumber in range(numberOfStructures):
        #   prior = single( priors( :, structureNumber ) ) / 65535;
        prior = priors[:, structureNumber] / 65535
        #   classNumber = reducingLookupTable( structureNumber );
        classNumber = reducingLookupTable[structureNumber]
        #
        #   likelihoods = zeros( length( maskIndices ), 1 );
        likelihoods = np.zeros(( likelihood_count, 1 ))
        #   numberOfComponents = numberOfGaussiansPerClass( classNumber );
        numberOfComponents = numberOfGaussiansPerClass[classNumber-1]
        #   for componentNumber = 1 : numberOfComponents
        for componentNumber in range(numberOfComponents):
            #     gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
            gaussianNumber = int(np.sum( numberOfGaussiansPerClass[: classNumber-1 ] ) + componentNumber)
            #
            #     mean = means( gaussianNumber, : )';v
            mean = ensure_dims(means, 2)[gaussianNumber, :].T
            #     variance = squeeze( variances( gaussianNumber, :, : ) );
            variance = ensure_dims(variances, 3)[gaussianNumber, :, : ]
            #     mixtureWeight = mixtureWeights( gaussianNumber );
            mixtureWeight = mixtureWeights[gaussianNumber]
            #
            #     L = chol( variance, 'lower' );  % variance = L * L'
            L = np.linalg.cholesky(variance)
            #     tmp = L \ ( data' - repmat( mean, [ 1 size( data, 1 ) ] ) );
            tmp = np.linalg.solve(L, data.T - mean)
            #     squaredMahalanobisDistances = ( sum( tmp.^2, 1 ) )';
            squaredMahalanobisDistances = ( np.sum( tmp**2, axis=0 ) ).T
            #     sqrtDeterminantOfVariance = prod( diag( L ) ); % Same as sqrt( det( variance ) )
            sqrtDeterminantOfVariance = np.prod( np.diag( L ) )
            #     gaussianLikelihoods = exp( -squaredMahalanobisDistances / 2 ) / ( 2 * pi )^( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance;
            gaussianLikelihoods = np.exp( -squaredMahalanobisDistances / 2 ) / ( 2 * np.pi )**( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance
            #
            #     likelihoods = likelihoods + gaussianLikelihoods * mixtureWeight;
            likelihoods = likelihoods + ensure_dims(gaussianLikelihoods, 2) * mixtureWeight
        #   end
        #
        #   posteriors( :, structureNumber ) = likelihoods .* prior;
        posteriors[:, structureNumber] = np.squeeze(likelihoods) * prior
        #
        # end % End loop over structures
    #
    # normalizer = sum( posteriors, 2 ) + eps;
    normalizer = np.sum( posteriors, 1 ) + eps
    # posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfStructures ] );
    posteriors = posteriors / ensure_dims(normalizer, 2)
    #
    #
    # % Compute volumes in mm^3
    # volumeOfOneVoxel = abs( det( imageToWorldTransformMatrix( 1:3, 1:3 ) ) );
    volumeOfOneVoxel = np.abs( np.linalg.det( imageToWorldTransformMatrix[0:3, 0:3] ) )
    # volumesInCubicMm = ( sum( posteriors ) )' * volumeOfOneVoxel;
    volumesInCubicMm = ( np.sum( posteriors, axis=0 ) ) * volumeOfOneVoxel
    #
    #
    # % Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
    # [ ~, structureNumbers ] = max( posteriors, [], 2 );
    structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)
    # freeSurferSegmentation = zeros( imageSize, 'uint16' );
    freeSurferSegmentation = np.zeros( imageSize, dtype=np.uint16 )
    # for structureNumber = 1 : numberOfStructures
        #   freeSurferSegmentation( maskIndices( find( structureNumbers == structureNumber ) ) ) = FreeSurferLabels( structureNumber );
    # end
    FreeSurferLabels = np.array(FreeSurferLabels, dtype=np.uint16)
    freeSurferSegmentation[mask == 1] = FreeSurferLabels[structureNumbers]
    # freeSurferSegmentation[mask == 1] = FreeSurferLabels[structureNumbers]

    # % Write to file, remembering to un-crop the segmentation to the original image size
    # uncroppedFreeSurferSegmentation = zeros( nonCroppedImageSize, 'single' );
    uncroppedFreeSurferSegmentation = np.zeros(nonCroppedImageSize, dtype=np.float32)
    # uncroppedFreeSurferSegmentation( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
    #                                  croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
    #                                  croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = freeSurferSegmentation;
    uncroppedFreeSurferSegmentation[ croppingOffset[0]: imageSize[0]+croppingOffset[0],
                                     croppingOffset[1]: imageSize[1]+croppingOffset[1],
                                     croppingOffset[2]: imageSize[2]+croppingOffset[2]] = freeSurferSegmentation
    # fprintf( 'Writing out freesurfer segmentation\n' );
    print('Writing out freesurfer segmentation')
    # kvlWriteImage( kvlCreateImage( uncroppedFreeSurferSegmentation ), ...
    #                fullfile( savePath, 'crispSegmentation.nii' ), ...
    #                imageToWorldTransform );
    GEMS2Python.KvlImage(require_np_array(uncroppedFreeSurferSegmentation)).write(
            os.path.join( savePath, 'crispSegmentation.nii' ),
            GEMS2Python.KvlTransform(require_np_array(imageToWorldTransformMatrix))
    )
    #
    #
    # % Also write out the bias field and the bias corrected image, each time remembering to un-crop the images
    # TODO: Are these bias field files even used other than debugging purposes?
    # for contrastNumber = 1 : numberOfContrasts
    #   [ dataPath, scanName, ext ] = fileparts( imageFileNames{ contrastNumber } );
    #
    #   % First bias field
    #   biasField = zeros( nonCroppedImageSize, 'single' );
    #   biasField( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
    #              croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
    #              croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = exp( biasFields( :, :, :, contrastNumber ) ) .* mask;
    #   outputFileName = fullfile( savePath, [ scanName '_biasField.nii' ] );
    #   kvlWriteImage( kvlCreateImage( biasField ), outputFileName, imageToWorldTransform );
    #
    #
    #   % Then bias field corrected data
    #   biasCorrected = zeros( nonCroppedImageSize, 'single' );
    #   biasCorrected( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
    #                  croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
    #                  croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = exp( biasCorrectedImageBuffers( :, :, :, contrastNumber ) );
    #   outputFileName = fullfile( savePath, [ scanName '_biasCorrected.nii' ] );
    #   kvlWriteImage( kvlCreateImage( biasCorrected ), outputFileName, imageToWorldTransform );
    #
    # end
    #
    #
    return {
        'FreeSurferLabels': FreeSurferLabels,
        'freeSurferSegmentation': freeSurferSegmentation,
        'uncroppedFreeSurferSegmentation': uncroppedFreeSurferSegmentation,
        'volumesInCubicMm': volumesInCubicMm,
    }

if __name__ == '__main__':
    checkpoint_manager = CheckpointManager()
    fixture = load_starting_fixture()
    part1_results_dict = checkpoint_manager.load('part1', 1)
    part2_results_dict = checkpoint_manager.load('part2', 1)
    part3_results_dict = samsegment_part3(
        fixture['modelSpecifications'],
        fixture['optimizationOptions'],
        part1_results_dict,
        part2_results_dict,
        checkpoint_manager
    )
    if checkpoint_manager:
        checkpoint_manager.save(part3_results_dict, 'part3', 1)
    names = part1_results_dict['names']
    FreeSurferLabels = part3_results_dict['FreeSurferLabels']
    volumesInCubicMm = part3_results_dict['volumesInCubicMm']
