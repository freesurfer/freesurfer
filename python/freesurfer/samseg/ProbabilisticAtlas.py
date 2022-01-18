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

        numberOfBlocks = 8
        if numberOfBlocks == 1:
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
            import time
            while True:
                globalTic = time.perf_counter()
                minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
                globalToc = time.perf_counter()
                print( f"  Total time spent: {globalToc-globalTic:0.4f} sec" )
                print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                historyOfMaximalDeformation.append(maximalDeformation)
                if maximalDeformation == 0:
                    break

        else:
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
                
                numberOfNodes = mesh.point_count
                numberOfNodesPerBlock = np.ceil( numberOfNodes / numberOfBlocks ).astype( 'int' )
                masks, submeshes, optimizers = [], [], []
                calculators = []
                import time
                useSmartClustering = True
                if useSmartClustering:
                    # Use k-means to group 
                    from sklearn.cluster import KMeans
                    nodePositions = mesh.points
                    kmeans = KMeans(n_clusters=numberOfBlocks, random_state=0).fit( nodePositions )
                for blockNumber in range( numberOfBlocks ):
                    tic = time.perf_counter()
                    start = blockNumber * numberOfNodesPerBlock
                    end = min( (blockNumber+1) * numberOfNodesPerBlock, numberOfNodes )
                    
                    if useSmartClustering:
                        mask = ( kmeans.labels_ == blockNumber )
                    else:  
                        mask = np.zeros( numberOfNodes, dtype=bool )
                        mask[ start:end ] = True
                    origMask = mask.copy()    
                    submesh = mesh.get_submesh( mask );
                    optimizer = gems.KvlOptimizer(optimizerType, submesh, calculator, optimizationParameters)
                    
                    masks.append( mask )
                    submeshes.append( submesh )
                    optimizers.append( optimizer )
                    toc = time.perf_counter()
                    
                    print( f"blockNumber {blockNumber}" )
                    print( f"     efficiency: {mask.sum() / origMask.sum() * 100:.4f} %" )
                    print( f"     relative size: {mask.sum() / numberOfNodes * 100:.4f} %" )
                    print( f"     setup time: {toc-tic:0.4f} sec" )
                  
                    
                
                if True:
                    # Save for reuse
                    self.optimizer = optimizers
                    self.submeshes = submeshes
                    self.masks = masks
                
            else:    
                print( "ProbabilisticAtlas reusing same optimizer!!" )
                optimizers = self.optimizer
                submeshes = self.submeshes
                masks = self.masks

                for optimizer in optimizers:
                    optimizer.set_calculator( calculator )
                

            # Run deformation optimization
            #historiesOfDeformationCost = [ [] for i in range( numberOfBlocks ) ]
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            
            #
            debug = False
            computeHistoryOfDeformationCost = False # Useful for analyzing convergence, but slow
            if computeHistoryOfDeformationCost: historyOfDeformationCost = []
            currentNodePositions = nodePositionsBeforeDeformation.copy()
            while True:
                import time
                globalTic = time.perf_counter()
                maximalDeformation = 0.0
                for blockNumber in range( numberOfBlocks ):
                    if debug:
                        mesh.points = currentNodePositions
                        tmpGlobalCostBefore, _ = calculator.evaluate_mesh_position( mesh )

                    submesh = submeshes[ blockNumber ]
                    mask = masks[ blockNumber ]
                    optimizer = optimizers[ blockNumber ]

                    submesh.points = currentNodePositions[ mask, : ]
                    
                    if debug:
                        tmpLocalCostBefore, _ = calculator.evaluate_mesh_position( submesh )
                    
                    blockMinLogLikelihoodTimesDeformationPrior, blockMaximalDeformation = optimizer.step_optimizer_samseg()
                    
                    currentNodePositions[ mask, : ] = submesh.points

                    if debug:
                        tmpLocalCostAfter, _ = calculator.evaluate_mesh_position( submesh )
                        mesh.points = currentNodePositions
                        tmpGlobalCostAfter, _ = calculator.evaluate_mesh_position( mesh )
                        print( f"  block {blockNumber}:" )
                        print( f"    local decrease: {tmpLocalCostBefore-tmpLocalCostAfter} ({tmpLocalCostBefore} - {tmpLocalCostAfter})" )
                        print( f"    global decrease: {tmpGlobalCostBefore-tmpGlobalCostAfter} ({tmpGlobalCostBefore} - {tmpGlobalCostAfter})" )

                    maximalDeformation = max( maximalDeformation, blockMaximalDeformation )
                    
                    # End loop over blocks
                    
                historyOfMaximalDeformation.append(maximalDeformation)

                globalToc = time.perf_counter()
                print( f"Total time spent: {globalToc-globalTic:0.4f} sec" )

                if computeHistoryOfDeformationCost:
                    tic = time.perf_counter()
                    mesh.points = currentNodePositions
                    minLogLikelihoodTimesDeformationPrior, _ = calculator.evaluate_mesh_position( mesh )
                    toc = time.perf_counter()
                    print( f'Additional time spent (unnecessarily) computing full cost function: {toc-tic:0.4f} sec' )
                    print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                    maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                    historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                else:
                    print( f"maximalDeformation={maximalDeformation:.4f}" )                  
                  
                    
                # Check early convergence    
                if maximalDeformation == 0:
                    break
                  
            #      
            mesh.points = currentNodePositions
            minLogLikelihoodTimesDeformationPrior, _ = calculator.evaluate_mesh_position( mesh )

            
            #
            if not computeHistoryOfDeformationCost: 
                historyOfDeformationCost = [ minLogLikelihoodTimesDeformationPrior ]

      
        # End test numberOfBlocks  

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
