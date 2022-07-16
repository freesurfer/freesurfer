import scipy
import numpy as np
from gems import KvlMeshCollection, KvlTransform, KvlOptimizer, KvlCostAndGradientCalculator, KvlImage


def kvlWarpMesh(sourceMeshCollectionFileName, sourceDeformation, targetMeshCollectionFileName):
    #
    # Applies sourceDeformation to the reference node positions of the sourceMeshCollection given by
    # sourceMeshCollectionFileName, and computes a targetDeformation on the reference node positions of the
    # targetMeshCollection given by targetMeshCollectionFileName that mimicks the same overall dense 3D deformation.
    #

    # Set flexibility of mesh deformation
    K = .001
    sourceMeshCollection = KvlMeshCollection()
    sourceMeshCollection.read(sourceMeshCollectionFileName)
    sourceReferenceMesh = sourceMeshCollection.reference_mesh
    sourceReferencePosition = sourceReferenceMesh.points
    sourceNumberOfNodes = sourceReferencePosition.shape[0]
    #identityTransform = KvlTransform(np.identity(4))
    targetMeshCollection = KvlMeshCollection()
    targetMeshCollection.read(targetMeshCollectionFileName)
    #targetMeshCollection.transform(identityTransform)
    targetMeshCollection.k = K
    targetReferenceMesh = targetMeshCollection.reference_mesh
    targetReferencePosition = targetReferenceMesh.points
    targetNumberOfNodes = targetReferencePosition.shape[0]
    
    # In the special case of identical mesh connectivity, no need to do any optimization
    if targetNumberOfNodes == sourceNumberOfNodes:
        if False:
            deltaReferencePosition = sourceReferencePosition - targetReferencePosition
            divergence = np.max(np.absolute(deltaReferencePosition))
            if divergence < 1e-2:
                # The reference meshes seem to be the same - therefore simply copy the deformation
                return sourceDeformation, 0.0, 0.0
        else:
            return sourceReferencePosition + sourceDeformation - targetReferencePosition, 0.0, 0.0
          
          
    imageSize = [int(1 + dim) for dim in np.max(sourceReferencePosition, axis=0)]
    # Rasterize the deformation by abusing alpha drawer
    deformation = np.copy(sourceDeformation)
    denseDeformation = np.zeros(imageSize + [3], dtype=np.double)
    if np.max(np.absolute(deformation)) > 0:
        maxDeformation = np.max(deformation)
        minDeformation = np.min(deformation)
        deltaDeformation = maxDeformation - minDeformation
        sourceReferenceMesh.alphas = (deformation - minDeformation) / deltaDeformation
        tmp = sourceReferenceMesh.rasterize_warp(imageSize, -1)
        # Unvisited voxels are marked by zeroes in all three coordinates - except possibly for the origin
        # which has all three coordinates zero as its natural state
        validMask = np.absolute(tmp)
        validMask = np.sum(validMask, axis=3)  ## only zero where all three coordinates are zero
        denseDeformation = (tmp) / (2 ** 16 - 1) * (maxDeformation - minDeformation) + minDeformation

        # Due to tetrahedral inside/outside checking performed in the rasterizer, some voxels on the boundary of the
        # image grid are never visited. There also seems to be non-visited voxels inside the image grid (bug in rasterizer).
        # We can resolve both issues by giving to each non-visited voxel the deformation field value of the closest voxel
        # that has been visited, augmented by knowledge that certain deformation components of voxels on the image grid
        # boundaries are zero by theory (sliding boundary conditions)
        closestIndices = scipy.ndimage.morphology.distance_transform_edt(
            validMask == 0, return_indices=True, return_distances=False)
        denseDeformation = denseDeformation[closestIndices[0], closestIndices[1], closestIndices[2]]

        # Boundaries are known a priori
        denseDeformation[0,:,:,0] = 0
        denseDeformation[-1,:,:,0] = 0
        denseDeformation[:,0,:,1] = 0
        denseDeformation[:,-1,:,1] = 0
        denseDeformation[:,:,0,2] = 0
        denseDeformation[:,:,-1,2] = 0

        ## Non zero values will have closestIndices pointing to themselves, so will not change
        ## Missed voxels with all 3 values at zero will map from nearest (euclidean) non-missed voxel
    # OK this seems to work. Now get (interpolate) the deformation at the target mesh nodes (in reference position)
    targetReferencePositionMap = np.transpose(targetReferencePosition)
    desiredTargetNodePosition = [
        scipy.ndimage.map_coordinates(denseDeformation[:, :, :, 0],
                                targetReferencePositionMap, mode='nearest', order=1, output=np.double),
        scipy.ndimage.map_coordinates(denseDeformation[:, :, :, 1],
                                targetReferencePositionMap, mode='nearest', order=1, output=np.double),
        scipy.ndimage.map_coordinates(denseDeformation[:, :, :, 2],
                                targetReferencePositionMap, mode='nearest', order=1, output=np.double),
    ]
    desiredTargetNodePosition += targetReferencePositionMap
    desiredTargetNodePosition = np.transpose(desiredTargetNodePosition)

    # Now deform the target mesh to try and bring the position of its mesh nodes close(r) to their target positions
    image = KvlImage(np.zeros((10, 10, 10), dtype=np.single))
    transform = KvlTransform(np.diag([1.0] * 4))
    calculator = KvlCostAndGradientCalculator('PointSet', [image], 'Sliding',
                                                          transform, [], [], [], [],
                                                          desiredTargetNodePosition
                                                          )
    # Get an optimizer, and stick the cost function into it
    optimizerType = 'L-BFGS'
    maximalDeformationStopCriterion = 0.005
    lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion
    optimizer = KvlOptimizer(optimizerType, targetReferenceMesh, calculator, {
        'Verbose': 1.0,
        'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
        'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
        'BFGS-MaximumMemoryLength': 12,
    })
    numberOfIterations = 0
    while True:
        minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_warp()
        if maximalDeformation == 0:
            break;
        numberOfIterations += 1
    targetNodePositions = targetReferenceMesh.points
    targetDeformation = targetNodePositions - targetReferencePosition
    targetDiscrepancy = targetNodePositions - desiredTargetNodePosition
    distances = np.sqrt(np.sum(targetDiscrepancy * targetDiscrepancy, axis=1))
    averageDistance = np.sum(distances) / distances.shape[0]
    maximumDistance = np.max(distances)
    return targetDeformation, averageDistance, maximumDistance
