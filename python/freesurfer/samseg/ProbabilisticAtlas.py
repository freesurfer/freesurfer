import numpy as np
from . import gems
from freesurfer.samseg.warp_mesh import kvlWarpMesh
from freesurfer.samseg.utilities import requireNumpyArray
import freesurfer as fs


class ProbabilisticAtlas:
    def __init__(self):
        # pass
        self.optimizer = None

    def getMesh(self, meshCollectionFileName,
                transform=None,
                K=None,
                initialDeformation=None, initialDeformationMeshCollectionFileName=None,
                returnInitialDeformationApplied=False):

        # Get the mesh
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(meshCollectionFileName)
        if K is not None:
            mesh_collection.k = K
        if transform:
            mesh_collection.transform(transform)
        else:
            transform = gems.KvlTransform(requireNumpyArray(np.eye(4)))
        mesh = mesh_collection.reference_mesh

        # See if we need to warp it
        estimatedNodeDeformationInTemplateSpace = None
        if initialDeformation is not None:
            #
            if initialDeformationMeshCollectionFileName is None:
                initialDeformationMeshCollectionFileName = meshCollectionFileName

            # Get the mesh node positions transformed back into template space (i.e., undoing the affine registration that we applied)
            nodePositions = mesh.points
            nodePositionsInTemplateSpace = self.mapPositionsFromSubjectToTemplateSpace(nodePositions, transform)

            # Get the estimated warp in template space
            [estimatedNodeDeformationInTemplateSpace, estimated_averageDistance,
             estimated_maximumDistance] = kvlWarpMesh(
                initialDeformationMeshCollectionFileName,
                initialDeformation,
                meshCollectionFileName
            )

            # Apply this warp on the mesh node positions in template space, and transform into current space
            desiredNodePositionsInTemplateSpace = nodePositionsInTemplateSpace + estimatedNodeDeformationInTemplateSpace
            desiredNodePositions = self.mapPositionsFromTemplateToSubjectSpace(desiredNodePositionsInTemplateSpace,
                                                                               transform)
            mesh.points = requireNumpyArray(desiredNodePositions)

        # Return what we got
        if returnInitialDeformationApplied:
            if estimatedNodeDeformationInTemplateSpace is None:
                estimatedNodeDeformationInTemplateSpace = np.zeros_like(mesh.points)
            return mesh, estimatedNodeDeformationInTemplateSpace
        else:
            return mesh

    def mapPositionsFromSubjectToTemplateSpace(self, positions, transform):

        #
        tmp = np.linalg.solve(transform.as_numpy_array,
                              np.pad(positions, ((0, 0), (0, 1)), mode='constant', constant_values=1).T).T
        return tmp[:, 0:3]

    def mapPositionsFromTemplateToSubjectSpace(self, positions, transform):

        #
        tmp = (transform.as_numpy_array @ \
               np.pad(positions, ((0, 0), (0, 1)), 'constant', constant_values=1).T).T
        return tmp[:, 0:3]

    def deformMesh(self, mesh, transform, data, mask, means, variances, mixtureWeights, numberOfGaussiansPerClass,
                   userOptimizationParameters={}):

        # Get images in ITK format
        numberOfContrasts = data.shape[-1]
        images = []
        for contrastNumber in range(numberOfContrasts):
            tmp = np.zeros(mask.shape, order='F')
            tmp[mask] = data[:, contrastNumber] + 1e-5 # Voxels with zero values but inside the mask 
                                                       # should not be skipped in the C++ code!
            images.append(gems.KvlImage(requireNumpyArray(tmp)))

        # Set up cost calculator
        calculator = gems.KvlCostAndGradientCalculator(
            typeName='AtlasMeshToIntensityImage',
            images=images,
            boundaryCondition='Sliding',
            transform=transform,
            means=means,
            variances=variances,
            mixtureWeights=mixtureWeights,
            numberOfGaussiansPerClass=numberOfGaussiansPerClass)

        # Get optimizer and plug calculator in it
        if self.optimizer is None:
            optimizerType = 'L-BFGS'
            #optimizerType = 'PartiallySeparable'
            optimizationParameters = {
                'Verbose': False,
                'MaximalDeformationStopCriterion': 0.001,  # measured in pixels,
                'LineSearchMaximalDeformationIntervalStopCriterion': 0.001,
                'MaximumNumberOfIterations': 20,
                'BFGS-MaximumMemoryLength': 12
            }
            optimizationParameters.update(userOptimizationParameters)
            optimizer = gems.KvlOptimizer(optimizerType, mesh, calculator, optimizationParameters)
            
            self.optimizer = optimizer
        else:    
            print( "ProbabilisticAtlas reusing same optimizer!!" )
            optimizer = self.optimizer
            optimizer.set_calculator( calculator )
            

        # Run deformation optimization
        historyOfDeformationCost = []
        historyOfMaximalDeformation = []
        nodePositionsBeforeDeformation = mesh.points
        
        #
        if True:
            numberOfNodes = mesh.point_count
            #mesh.can_moves.sum(axis=0) / numberOfNodes                                                                                          
            orig_can_moves = mesh.can_moves                                                                                                     
            numberOfBlocks = 2                                                                                                                  
            numberOfNodesPerBlock = np.ceil( numberOfNodes / numberOfBlocks ).astype( 'int' )                                                                  
            blockNumber=1                                                                                                                      
            start = blockNumber * numberOfNodesPerBlock
            end = min( (blockNumber+1) * numberOfNodesPerBlock, numberOfNodes )
            can_moves = np.zeros_like( orig_can_moves )
            can_moves[ start:end, : ] = orig_can_moves[ start:end, : ]
            mesh.can_moves = can_moves
        
        
        while True:
            minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
            print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
            maximalDeformation, minLogLikelihoodTimesDeformationPrior))
            historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
            historyOfMaximalDeformation.append(maximalDeformation)
            if maximalDeformation == 0:
                break

        # Return
        nodePositionsAfterDeformation = mesh.points
        maximalDeformationApplied = np.sqrt(
            np.max(np.sum((nodePositionsAfterDeformation - nodePositionsBeforeDeformation) ** 2, 1)))
        return historyOfDeformationCost, historyOfMaximalDeformation, maximalDeformationApplied, minLogLikelihoodTimesDeformationPrior

    def saveDeformedAtlas(self, originalAtlasFileName, deformedAtlasFileName, arg, applyAsDeformation=False):

        #
        self.mesh_collection = gems.KvlMeshCollection()
        self.mesh_collection.read(originalAtlasFileName)
        if not applyAsDeformation:
            position = arg
        else:
            position = self.mesh_collection.reference_mesh.points + arg
        self.mesh_collection.reference_mesh.points = requireNumpyArray(position)
        self.mesh_collection.write(deformedAtlasFileName)
