import traceback

from scipy import ndimage

import GEMS2Python
import numpy as np
import scipy.ndimage
import scipy.io


# function [ targetDeformation, averageDistance, maximumDistance ] = kvlWarpMesh( sourceMeshCollectionFileName, sourceDeformation, targetMeshCollectionFileName, showFigures )
def kvlWarpMesh(sourceMeshCollectionFileName, sourceDeformation, targetMeshCollectionFileName, showFigures=False):
    # %
    # % Applies sourceDeformation to the reference node positions of the sourceMeshCollection given by
    # % sourceMeshCollectionFileName, and computes a targetDeformation on the reference node positions of the
    # % targetMeshCollection given by targetMeshCollectionFileName that mimicks the same overall dense 3D deformation.
    # %
    #
    #
    # if ( nargin == 0 )
    #   %
    #   % Test run this function
    #   %
    #
    #   %
    #   if 1
    #     sourceMeshCollectionFileName = ...
    #         '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30_multires6.gz';
    #   else
    #     % Test special case of matching mesh connectivity
    #     sourceMeshCollectionFileName = ...
    #         '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30.gz';
    #   end
    #   targetMeshCollectionFileName = ...
    #       '/data/tmp/tmpBuckner/tmp3/scratch/koenLogDir10SubjectAtlas3/CurrentMeshCollection30.gz';
    #
    #
    #   %
    #   sourceMeshCollection = kvlReadMeshCollection( sourceMeshCollectionFileName );
    #   sourceReferenceMesh = kvlGetMesh( sourceMeshCollection, -1 );
    #   sourceReferencePosition = kvlGetMeshNodePositions( sourceReferenceMesh );
    #   sourceDeformedMesh = kvlGetMesh( sourceMeshCollection, 8 ); % Purposefully deformed
    #   sourceDeformedPosition = kvlGetMeshNodePositions( sourceDeformedMesh );
    #   sourceDeformation = sourceDeformedPosition - sourceReferencePosition;
    #
    #   %
    #   showFigures = true;
    #   targetDeformation = kvlWarpMesh( sourceMeshCollectionFileName, sourceDeformation, ...
    #                                    targetMeshCollectionFileName, showFigures );
    #
    #   % Write deformed target mesh collection to file
    #   targetMeshCollection = kvlReadMeshCollection( targetMeshCollectionFileName );
    #   targetReferenceMesh = kvlGetMesh( targetMeshCollection, -1 );
    #   targetReferencePosition = kvlGetMeshNodePositions( targetReferenceMesh );
    #   targetPosition = targetReferencePosition + targetDeformation;
    #   kvlSetMeshCollectionPositions( targetMeshCollection, targetReferencePosition, targetPosition );
    #   [a b c ] = fileparts( targetMeshCollectionFileName );
    #   resultFileName = fullfile( a, [ b '_debugDeformed' c ] );
    #   kvlWriteMeshCollection( targetMeshCollection, resultFileName );
    #   disp( [ 'Wrote result in file: ' resultFileName ] )
    #
    #   return;
    # end
    #
    # if ( nargin < 4 )
    #   showFigures = false;
    # end
    #
    #
    # % Set flexibility of mesh deformation
    # K = .001;
    K = .001
    #
    #
    # %
    # sourceMeshCollection = kvlReadMeshCollection( sourceMeshCollectionFileName );
    sourceMeshCollection = GEMS2Python.KvlMeshCollection()
    sourceMeshCollection.read(sourceMeshCollectionFileName)
    # sourceReferenceMesh = kvlGetMesh( sourceMeshCollection, -1 );
    sourceReferenceMesh = sourceMeshCollection.reference_mesh
    # sourceReferencePosition = kvlGetMeshNodePositions( sourceReferenceMesh );
    sourceReferencePosition = sourceReferenceMesh.points
    # sourceNumberOfNodes = size( sourceReferencePosition, 1 );
    sourceNumberOfNodes = sourceReferencePosition.shape[0]
    #
    # %
    identityTransform = GEMS2Python.KvlTransform(np.identity(4))
    # targetMeshCollection = kvlReadMeshCollection( targetMeshCollectionFileName, kvlCreateTransform( eye( 4 ) ), K );
    targetMeshCollection = GEMS2Python.KvlMeshCollection()
    targetMeshCollection.read(targetMeshCollectionFileName)
    targetMeshCollection.transform(identityTransform)
    targetMeshCollection.k = K
    # targetReferenceMesh = kvlGetMesh( targetReferenceMesh, -1 );
    targetReferenceMesh = targetMeshCollection.reference_mesh
    # targetReferencePosition = kvlGetMeshNodePositions( targetReferenceMesh );
    targetReferencePosition = targetReferenceMesh.points
    # targetNumberOfNodes = size( targetReferencePosition, 1 );
    targetNumberOfNodes = targetReferencePosition.shape[0]
    #
    # % In the special case of identical mesh connectivity, no need to do any optimization
    # if ( targetNumberOfNodes == sourceNumberOfNodes )
    if targetNumberOfNodes == sourceNumberOfNodes:
        #   if ( max( abs( sourceReferencePosition(:) - targetReferencePosition(:) ) ) < 1e-2 )
        deltaReferencePosition = sourceReferencePosition - targetReferencePosition
        divergence = np.max(np.absolute(deltaReferencePosition))
        if divergence < 1e-2:
            #     % The reference meshes seem to be the same - therefore simply copy the deformation
            #     targetDeformation = sourceDeformation;
            #     averageDistance = 0.0;
            #     maximumDistance = 0.0;
            #
            #     return
            return sourceDeformation, 0.0, 0.0
    #   end
    # end
    #
    #
    # %
    # imageSize = max( sourceReferencePosition ) + 1;
    imageSize = [int(1 + dim) for dim in np.max(sourceReferencePosition, axis=0)]
    #
    #
    # % Rasterize the deformation by abusing alpha drawer
    # deformation = sourceDeformation;
    deformation = np.copy(sourceDeformation)
    # denseDeformation = zeros( [ imageSize 3 ] );
    denseDeformation = np.zeros(imageSize + [3], dtype=np.double)
    # if ( max( abs( deformation(:) ) ) > 0 )
    if np.max(np.absolute(deformation)) > 0:
        #   %
        #   maxDeformation = max( deformation(:) );
        maxDeformation = np.max(deformation, axis=0)
        #   minDeformation = min( deformation(:) );
        minDeformation = np.min(deformation, axis=0)
        #
        #   kvlSetAlphasInMeshNodes( sourceReferenceMesh, single( deformation - minDeformation ) ./ ...
        #                                                 repmat( maxDeformation - minDeformation, [ sourceNumberOfNodes 3 ] ) );
        deltaDeformation = maxDeformation - minDeformation
        sourceReferenceMesh.alphas = (deformation - minDeformation) / deltaDeformation
        #   % denseDeformation = double( kvlRasterizeAtlasMesh( sourceReferenceMesh, imageSize ) ) / ( 2^16 -1 ) ...
        #   %                               * ( maxDeformation - minDeformation ) + minDeformation;
        #   tmp = kvlRasterizeAtlasMesh( sourceReferenceMesh, imageSize );
        tmp = sourceReferenceMesh.rasterize(imageSize, -1)
        #
        #   % Unvisited voxels are marked by zeroes in all three coordinates - except possibly for the origin
        #   % which has all three coordinates zero as its natural state
        #   validMask = ( sum( tmp ~= 0, 4 ) ~= 0 );
        validMask = np.absolute(tmp)
        validMask = np.sum(validMask, axis=3)  ## only zero where all three coordinates are zero
        #   validMask( 1, 1, 1 ) = 1;
        validMask[0, 0, 0] = 1  ## origin is allowed to be zero
        #   if showFigures
        #     figure
        #     showImage( validMask );
        #   end
        #
        #   %
        #   denseDeformation = double( tmp ) / ( 2^16 -1 ) * ( maxDeformation - minDeformation ) + minDeformation;
        denseDeformation = (tmp) / (2 ** 16 - 1) * (maxDeformation - minDeformation) + minDeformation
        #
        #   % Due to tetrahedral inside/outside checking performed in the rasterizer, some voxels on the boundary of the
        #   % image grid are never visited. There also seems to be non-visited voxels inside the image grid (bug in rasterizer).
        #   % We can resolve both issues by giving to each non-visited voxel the deformation field value of the closest voxel
        #   % that has been visited, augmented by knowledge that certain deformation components of voxels on the image grid
        #   % boundaries are zero by theory (sliding boundary conditions)
        #
        #   [ distances, closestIndices ] = bwdist( validMask );
        closestIndices = scipy.ndimage.morphology.distance_transform_edt(
            validMask == 0, return_indices=True, return_distances=False)
        #   for directionNumber = 1 : 3
        #     corrected = denseDeformation( :, :, :, directionNumber );
        #     corrected( find( ~validMask ) ) =  corrected( closestIndices( ~validMask ) );
        #
        #     % Boundaries are known a priori
        #     if ( directionNumber == 1 )
        #       corrected( 1, :, : ) = 0;
        #       corrected( end, :, : ) = 0;
        #     elseif ( directionNumber == 2 )
        #       corrected( :, 1, : ) = 0;
        #       corrected( :, end, : ) = 0;
        #     else
        #       corrected( :, :, 1 ) = 0;
        #       corrected( :, :, end ) = 0;
        #     end
        #
        #     denseDeformation( :, :, :, directionNumber ) = corrected;
        #
        #   end % End loop over directions
        #
        # end % End test if there is deformation at all
        ## Non zero values will have closestIndices pointing to themselves, so will not change
        ## Missed voxels with all 3 values at zero will map from nearest (euclidean) non-missed voxel
        # TODO: make this efficient
        # newDenseDeformation = denseDeformation[closestIndices]
        if False:
            x_map = closestIndices[0]
            y_map = closestIndices[1]
            z_map = closestIndices[2]
            x_limit, y_limit, z_limit, should_be_three = denseDeformation.shape
            for x in range(x_limit):
                for y in range(y_limit):
                    for z in range(z_limit):
                        close_x = x_map[x, y, z]
                        close_y = y_map[x, y, z]
                        close_z = z_map[x, y, z]
                        if x != close_x or y != close_y or z != close_z:
                            denseDeformation[x, y, z] = denseDeformation[close_x, close_y, close_z]
    #
    # %
    # if showFigures
    #   for dimensionNumber = 1 : 3
    #     figure
    #     showImage( denseDeformation( :, :, :, dimensionNumber ) )
    #     drawnow
    #   end
    # end
    #
    #
    # %
    # if showFigures
    #   [ x y z ] = ndgrid( [ 0 : imageSize(1)-1 ], [ 0 : imageSize(2)-1 ], [ 0 : imageSize(3)-1 ] );
    #   densePositions = zeros( [ imageSize 3 ] );
    #   densePositions( :, :, :, 1 ) = x;
    #   densePositions( :, :, :, 2 ) = y;
    #   densePositions( :, :, :, 3 ) = z;
    #
    #   figure
    #   q = quiver3( densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 1 ), ...
    #                densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 2 ), ...
    #                densePositions( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 3 ), ...
    #                denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 1 ), ...
    #                denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 2 ), ...
    #                denseDeformation( 1 : 20 : end, 1 : 20 : end, 1 : 20 : end, 3 ) );
    #    title( 'Dense deformation' )
    #
    # end
    #
    # % OK this seems to work. Now get (interpolate) the deformation at the target mesh nodes (in reference position)
    # desiredTargetNodePosition = zeros( targetNumberOfNodes, 3 );
    # for directionNumber = 1 : 3
    #   desiredTargetNodePosition( :, directionNumber ) = ...
    #             targetReferencePosition( :, directionNumber ) + ...
    #             interpn( denseDeformation( :, :, :, directionNumber ), ...
    #                      targetReferencePosition( :, 1 )+1, ...
    #                      targetReferencePosition( :, 2 )+1, ...
    #                      targetReferencePosition( :, 3 )+1 );
    # end
    targetReferencePositionMap = np.transpose(targetReferencePosition)
    desiredTargetNodePosition = [
        ndimage.map_coordinates(denseDeformation[:, :, :, 0],
                                targetReferencePositionMap, mode='nearest', output=np.double),
        ndimage.map_coordinates(denseDeformation[:, :, :, 1],
                                targetReferencePositionMap, mode='nearest', output=np.double),
        ndimage.map_coordinates(denseDeformation[:, :, :, 2],
                                targetReferencePositionMap, mode='nearest', output=np.double),
    ]
    desiredTargetNodePosition = np.transpose(desiredTargetNodePosition)
    #
    # %
    # if showFigures
    #   figure
    #   desiredDeformation = desiredTargetNodePosition - targetReferencePosition;
    #   quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
    #           targetReferencePosition( 1 : 10 : end, 2 ), ...
    #           targetReferencePosition( 1 : 10 : end, 3 ), ...
    #           desiredDeformation( 1 : 10 : end, 1 ), ...
    #           desiredDeformation( 1 : 10 : end, 2 ), ...
    #           desiredDeformation( 1 : 10 : end, 3 ) );
    #   title( 'desired deformation' )
    # end
    #
    #
    # % Now deform the target mesh to try and bring the position of its mesh nodes close(r) to their target positions
    # image = kvlCreateImage( zeros( 10, 10, 10, 'single' ) ); % Dummy but we need it with the current interface
    image = GEMS2Python.KvlImage(np.zeros((10, 10, 10), dtype=np.single))
    # transform = kvlCreateTransform( eye( 4 ) );
    transform = GEMS2Python.KvlTransform(np.diag([1.0] * 4))
    # calculator = kvlGetCostAndGradientCalculator( 'PointSet', image, 'Sliding', transform, [], [], [], [], ...
    #                                               desiredTargetNodePosition );
    calculator = GEMS2Python.KvlCostAndGradientCalculator('PointSet', [image], 'Sliding',
                                                          ## TODO: need full calculator
                                                          transform, [], [], [], [],
                                                          desiredTargetNodePosition
                                                          )
    #
    # % Get an optimizer, and stick the cost function into it
    # optimizerType = 'L-BFGS'; % 'FixedStepGradientDescent','GradientDescent','ConjugateGradient', or 'L-BFGS'
    optimizerType = 'L-BFGS'
    # maximalDeformationStopCriterion = 0.005; % Measured in voxels
    maximalDeformationStopCriterion = 0.005
    # lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion; % Doesn't seem to matter very much
    lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion
    # optimizer = kvlGetOptimizer( optimizerType, targetReferenceMesh, calculator, ...
    #                                 'Verbose', 1, ...
    #                                 'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
    #                                 'LineSearchMaximalDeformationIntervalStopCriterion', ...
    #                                 lineSearchMaximalDeformationIntervalStopCriterion, ...
    #                                 'BFGS-MaximumMemoryLength', 12 ); % Affine registration only has 12 DOF
    optimizer = GEMS2Python.KvlOptimizer(optimizerType, targetReferenceMesh, calculator, {
        'Verbose': 1.0,
        'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
        'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
        'BFGS-MaximumMemoryLength': 12,
    })
    #
    # numberOfIterations = 0;
    numberOfIterations = 0
    # startTime = tic;
    # if showFigures
    #   historyOfAverageDistance = [];
    #   historyOfMaximumDistance = [];
    #   iterationFigure = figure;
    # end
    # while true
    while True:
        #
        #   %
        #   if showFigures
        #     targetNodePositions = kvlGetMeshNodePositions( targetReferenceMesh );
        #     distances = sqrt( sum( ( targetNodePositions - desiredTargetNodePosition ).^2, 2 ) );
        #     averageDistance = sum( distances ) / length( distances )
        #     maximumDistance = max( distances );
        #     historyOfAverageDistance = [ historyOfAverageDistance; averageDistance ];
        #     historyOfMaximumDistance = [ historyOfMaximumDistance; maximumDistance ];
        #
        #     figure( iterationFigure )
        #     subplot( 2, 1, 1 )
        #     plot( historyOfAverageDistance )
        #     grid
        #     title( [ 'average distance: ' num2str( historyOfAverageDistance( end ) ) ] )
        #     subplot( 2, 1, 2 )
        #     plot( historyOfMaximumDistance )
        #     grid
        #     title( [ 'max distance: ' num2str( historyOfMaximumDistance( end ) ) ] )
        #     drawnow
        #   end
        #
        #   %
        #   [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
        minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer()
        #   %return
        #   if ( maximalDeformation == 0 )
        #     break;
        #   end
        if maximalDeformation == 0:
            break;
        #   numberOfIterations = numberOfIterations + 1;
        numberOfIterations += 1
        #
        # end
    # toc( startTime )
    #
    # %
    # targetNodePositions = kvlGetMeshNodePositions( targetReferenceMesh );
    targetNodePositions = targetReferenceMesh.points
    # targetDeformation = targetNodePositions - targetReferencePosition;
    targetDeformation = targetNodePositions - targetReferencePosition
    # distances = sqrt( sum( ( targetNodePositions - desiredTargetNodePosition ).^2, 2 ) );
    targetDiscrepancy = targetNodePositions - desiredTargetNodePosition
    distances = np.sqrt(np.sum(targetDiscrepancy * targetDiscrepancy, axis=1))
    # averageDistance = sum( distances ) / length( distances );
    averageDistance = np.sum(distances) / distances.shape[0]
    # maximumDistance = max( distances );
    maximumDistance = np.max(distances)
    #
    # if showFigures
    #   figure
    #   quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
    #            targetReferencePosition( 1 : 10 : end, 2 ), ...
    #            targetReferencePosition( 1 : 10 : end, 3 ), ...
    #            targetDeformation( 1 : 10 : end, 1 ), ...
    #            targetDeformation( 1 : 10 : end, 2 ), ...
    #            targetDeformation( 1 : 10 : end, 3 ) );
    #   title( 'Achieved deformation' )
    #
    #   figure
    #   nonAchievedDeformationComponent = targetDeformation - desiredDeformation;
    #   quiver3( targetReferencePosition( 1 : 10 : end, 1 ), ...
    #            targetReferencePosition( 1 : 10 : end, 2 ), ...
    #            targetReferencePosition( 1 : 10 : end, 3 ), ...
    #            nonAchievedDeformationComponent( 1 : 10 : end, 1 ), ...
    #            nonAchievedDeformationComponent( 1 : 10 : end, 2 ), ...
    #            nonAchievedDeformationComponent( 1 : 10 : end, 3 ), ...
    #            0 ); % Don't scale arrows
    #   title( 'Difference between achieved and desired deformation (true scale)' )
    #
    #
    # end
    #
    #
    #
    return targetDeformation, averageDistance, maximumDistance


if __name__ == '__main__':
    try:
        TEST_PATH = '/home/willy/work/cm/my_tests/kvlWarpMesh/kvlWarpMesh_'
        INPUT_PATH = TEST_PATH + 'input.mat'
        OUTPUT_PATH = TEST_PATH + 'output.mat'
        input_fixture = scipy.io.loadmat(INPUT_PATH, struct_as_record=False, squeeze_me=True)
        output_fixture = scipy.io.loadmat(OUTPUT_PATH, struct_as_record=False, squeeze_me=True)
        sourceDeformation = input_fixture['sourceDeformation']
        sourceMeshCollectionFileName = '/home/willy/work/cm/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlas_level1.txt.gz'
        targetMeshCollectionFileName = '/home/willy/work/cm/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlas_level2.txt.gz'
        targetDeformation, averageDistance, maximumDistance = kvlWarpMesh(sourceMeshCollectionFileName, sourceDeformation, targetMeshCollectionFileName)
        print("averageDistance = {0}".format(averageDistance))
        print("maximumDistance = {0}".format(maximumDistance))
        print("targetDeformation.shape = {0}".format(targetDeformation.shape))
    except Exception as flaw:
        print('flaw = {0}'.format(str(flaw)))
        traceback.print_exc()
