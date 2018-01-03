import scipy.io
import numpy as np
import os
import GEMS2Python
from dotdict import DotDict

from as_python.samseg.dev_utils.debug_client import request_var, compare_vars, CheckpointManager
from as_python.samseg.kvl_merge_alphas import kvlMergeAlphas
from as_python.samseg.bias_correction import backprojectKroneckerProductBasisFunctions

MATLAB_FIXTURE_PATH = '/Users/ys/work/freesurfer/GEMS2/Testing/matlab_data/'

fixture = scipy.io.loadmat(os.path.join(MATLAB_FIXTURE_PATH, 'part3.mat'), struct_as_record=False, squeeze_me=True)
locals().update(fixture)

eps = np.finfo(float).eps

checkpoint_manager = CheckpointManager()

def ensure_dims(np_array, dims):
    if np_array.ndim < dims:
        return np.expand_dims(np_array, axis=dims)
    elif np_array.ndim == dims:
        return np_array

def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])

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

transform = GEMS2Python.KvlTransform(require_np_array(transformMatrix))
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
# TODO: Need to remove this fixture once kvlWarpMesh
# estimatedNodeDeformationInTemplateSpace = ...
#           kvlWarpMesh( optimizationOptions.multiResolutionSpecification( end ).atlasFileName, ...
#                         historyWithinEachMultiResolutionLevel( end ).finalNodePositionsInTemplateSpace ...
#                         - historyWithinEachMultiResolutionLevel( end ).initialNodePositionsInTemplateSpace, ...
#                         modelSpecifications.atlasFileName );
#
# % Apply this warp on the mesh node positions in template space, and transform into current space
estimatedNodeDeformationInTemplateSpace = checkpoint_manager.load('estimatedNodeDeformationInTemplateSpace', 1)['estimatedNodeDeformationInTemplateSpace']
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
data = biasCorrectedImageBuffers.reshape( [np.prod( imageSize ), numberOfContrasts ], order='F' )
# priors = kvlRasterizeAtlasMesh( mesh, imageSize );
priors = mesh.rasterize(imageSize, -1)
# priors = reshape( priors, [ prod( imageSize ) numberOfStructures ] );
priors = priors.reshape([ np.prod( imageSize ), numberOfStructures ] , order='F' )
#
#
# % Ignore everything that's has zero intensity
# priors = priors( maskIndices, : );
priors = priors[maskIndices, :]
# data = data( maskIndices, : );
data = data[maskIndices, :]
#
#
# % Calculate the posteriors
# posteriors = zeros( size( priors ), 'double' );
posteriors = np.zeros_like(priors)
# for structureNumber = 1 : numberOfStructures
for structureNumber in range(numberOfStructures):
    #   prior = single( priors( :, structureNumber ) ) / 65535;
    prior = priors[:, structureNumber] / 65535
    #   classNumber = reducingLookupTable( structureNumber );
    classNumber = reducingLookupTable[structureNumber]
    #
    #   likelihoods = zeros( length( maskIndices ), 1 );
    likelihoods = np.zeros(( len( maskIndices ), 1 ))
    #   numberOfComponents = numberOfGaussiansPerClass( classNumber );
    numberOfComponents = numberOfGaussiansPerClass[classNumber]
    #   for componentNumber = 1 : numberOfComponents
    for componentNumber in range(numberOfComponents):
        #     gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
        gaussianNumber = np.sum( numberOfGaussiansPerClass[: classNumber-1 ] ) + componentNumber
        #
        #     mean = means( gaussianNumber, : )';v
        mean = means[gaussianNumber, : ].T
        #     variance = squeeze( variances( gaussianNumber, :, : ) );
        variance = np.squeeze( variances[gaussianNumber, :, : ] )
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
        likelihoods = likelihoods + gaussianLikelihoods * mixtureWeight
    #   end
    #
    #   posteriors( :, structureNumber ) = likelihoods .* prior;
    posteriors[:, structureNumber] = likelihoods * prior
    #
    # end % End loop over structures
#
# normalizer = sum( posteriors, 2 ) + eps;
normalizer = np.sum( posteriors, 2 ) + eps
# posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfStructures ] );
posteriors = posteriors / normalizer
#
#
# % Compute volumes in mm^3
# volumeOfOneVoxel = abs( det( imageToWorldTransformMatrix( 1:3, 1:3 ) ) );
volumeOfOneVoxel = np.abs( np.det( imageToWorldTransformMatrix[0:3, 0:3] ) )
# volumesInCubicMm = ( sum( posteriors ) )' * volumeOfOneVoxel;
volumesInCubicMm = ( np.sum( posteriors ) ).T @ volumeOfOneVoxel
#
#
# % Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
# [ ~, structureNumbers ] = max( posteriors, [], 2 );
structureNumbers = np.max( posteriors, [], 2 )
# freeSurferSegmentation = zeros( imageSize, 'uint16' );
freeSurferSegmentation = np.zeros( imageSize, dtype=np.uint16 )
# for structureNumber = 1 : numberOfStructures
for structureNumber in range(numberOfStructures):
    freeSurferSegmentation[maskIndices[structureNumbers == structureNumber]] = FreeSurferLabels( structureNumber )
    #   freeSurferSegmentation( maskIndices( find( structureNumbers == structureNumber ) ) ) = FreeSurferLabels( structureNumber );
# end
#
# % Write to file, remembering to un-crop the segmentation to the original image size
# uncroppedFreeSurferSegmentation = zeros( nonCroppedImageSize, 'single' );
# uncroppedFreeSurferSegmentation( croppingOffset( 1 ) + [ 1 : imageSize( 1 ) ], ...
#                                  croppingOffset( 2 ) + [ 1 : imageSize( 2 ) ], ...
#                                  croppingOffset( 3 ) + [ 1 : imageSize( 3 ) ] ) = freeSurferSegmentation;
# fprintf( 'Writing out freesurfer segmentation\n' );
# kvlWriteImage( kvlCreateImage( uncroppedFreeSurferSegmentation ), ...
#                fullfile( savePath, 'crispSegmentation.nii' ), ...
#                imageToWorldTransform );
#
#
# % Also write out the bias field and the bias corrected image, each time remembering to un-crop the images
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
