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
# biasFields = zeros( [ imageSize numberOfContrasts ] );
# for contrastNumber = 1 : numberOfContrasts
#   biasField = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, biasFieldCoefficients( :, contrastNumber ) );
#   biasCorrectedImageBuffers( :, :, :, contrastNumber ) = ...
#               imageBuffers( :, :, :, contrastNumber ) - biasField .* mask;
#   biasFields( :, :, :, contrastNumber ) = biasField;
# end
#
#
#
# % Read the atlas, applying the affine registration transform
# meshCollection = kvlReadMeshCollection( modelSpecifications.atlasFileName, transform, modelSpecifications.K );
# mesh = kvlGetMesh( meshCollection, -1 );
#
# % Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
# nodePositions = kvlGetMeshNodePositions( mesh );
# numberOfNodes = size( nodePositions, 1 );
# transformMatrix = double( kvlGetTransformMatrix( transform ) );
# tmp = ( transformMatrix \ [ nodePositions ones( numberOfNodes, 1 ) ]' )';
# nodePositionsInTemplateSpace = tmp( :, 1 : 3 );
#
# % Get the estimated warp in template space
# estimatedNodeDeformationInTemplateSpace = ...
#           kvlWarpMesh( optimizationOptions.multiResolutionSpecification( end ).atlasFileName, ...
#                         historyWithinEachMultiResolutionLevel( end ).finalNodePositionsInTemplateSpace ...
#                         - historyWithinEachMultiResolutionLevel( end ).initialNodePositionsInTemplateSpace, ...
#                         modelSpecifications.atlasFileName );
#
# % Apply this warp on the mesh node positions in template space, and transform into current space
# desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace;
# tmp = ( transformMatrix * [ desiredNodePositionsInTemplateSpace ones( numberOfNodes, 1 ) ]' )';
# desiredNodePositions = tmp( :, 1 : 3 );
#
# %
# kvlSetMeshNodePositions( mesh, desiredNodePositions );
#
# %
# alphas = kvlGetAlphasInMeshNodes( mesh );
# numberOfStructures = size( alphas, 2 );
#
#
#
# % Get the priors as dictated by the current mesh position
# data = reshape( biasCorrectedImageBuffers, [ prod( imageSize ) numberOfContrasts ] );
# priors = kvlRasterizeAtlasMesh( mesh, imageSize );
# priors = reshape( priors, [ prod( imageSize ) numberOfStructures ] );
#
#
# % Ignore everything that's has zero intensity
# priors = priors( maskIndices, : );
# data = data( maskIndices, : );
#
#
# % Calculate the posteriors
# posteriors = zeros( size( priors ), 'double' );
# for structureNumber = 1 : numberOfStructures
#
#   prior = single( priors( :, structureNumber ) ) / 65535;
#   classNumber = reducingLookupTable( structureNumber );
#
#   likelihoods = zeros( length( maskIndices ), 1 );
#   numberOfComponents = numberOfGaussiansPerClass( classNumber );
#   for componentNumber = 1 : numberOfComponents
#     gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
#
#     mean = means( gaussianNumber, : )';v
#     variance = squeeze( variances( gaussianNumber, :, : ) );
#     mixtureWeight = mixtureWeights( gaussianNumber );
#
#     L = chol( variance, 'lower' );  % variance = L * L'
#     tmp = L \ ( data' - repmat( mean, [ 1 size( data, 1 ) ] ) );
#     squaredMahalanobisDistances = ( sum( tmp.^2, 1 ) )';
#     sqrtDeterminantOfVariance = prod( diag( L ) ); % Same as sqrt( det( variance ) )
#     gaussianLikelihoods = exp( -squaredMahalanobisDistances / 2 ) / ( 2 * pi )^( numberOfContrasts / 2 ) / sqrtDeterminantOfVariance;
#
#     likelihoods = likelihoods + gaussianLikelihoods * mixtureWeight;
#   end
#
#   posteriors( :, structureNumber ) = likelihoods .* prior;
#
# end % End loop over structures
#
# normalizer = sum( posteriors, 2 ) + eps;
# posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfStructures ] );
#
#
# % Compute volumes in mm^3
# volumeOfOneVoxel = abs( det( imageToWorldTransformMatrix( 1:3, 1:3 ) ) );
# volumesInCubicMm = ( sum( posteriors ) )' * volumeOfOneVoxel;
#
#
# % Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
# [ ~, structureNumbers ] = max( posteriors, [], 2 );
# freeSurferSegmentation = zeros( imageSize, 'uint16' );
# for structureNumber = 1 : numberOfStructures
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
