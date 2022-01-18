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

        numberOfBlocks = 4
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

        elif False:
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
                
                optimizers = []
                for blockNumber in range( 0, numberOfBlocks ): 
                    optimizer = gems.KvlOptimizer(optimizerType, mesh, calculator, optimizationParameters)
                    optimizers.append( optimizer )
                
                self.optimizer = optimizers
            else:    
                print( "ProbabilisticAtlas reusing same optimizer!!" )
                optimizers = self.optimizer
                for optimizer in optimizers:
                    optimizer.set_calculator( calculator )
                

            # Run deformation optimization
            historyOfDeformationCost = []
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            
            #
            numberOfNodes = mesh.point_count
            #mesh.can_moves.sum(axis=0) / numberOfNodes                                                                                          
            orig_can_moves = mesh.can_moves                                                                                                     
            numberOfNodesPerBlock = np.ceil( numberOfNodes / numberOfBlocks ).astype( 'int' )
            while True:
                #timeSpentDoingCxxCall = 0;
                import time
                globalTic = time.perf_counter()
                maximalDeformation = 0.0
                for blockNumber in range( numberOfBlocks ):
                    if False:
                        mesh.can_moves = orig_can_moves
                        tmpGlobalCostBefore, _ = calculator.evaluate_mesh_position( mesh )
                    start = blockNumber * numberOfNodesPerBlock
                    end = min( (blockNumber+1) * numberOfNodesPerBlock, numberOfNodes )
                    can_moves = np.zeros_like( orig_can_moves )
                    can_moves[ start:end, : ] = orig_can_moves[ start:end, : ]
                    mesh.can_moves = can_moves
                    if False:
                        tmpLocalCostBefore, _ = calculator.evaluate_mesh_position( mesh )
                    optimizer = optimizers[ blockNumber ]
                    #tic = time.perf_counter()
                    blockMinLogLikelihoodTimesDeformationPrior, blockMaximalDeformation = optimizer.step_optimizer_samseg()
                    #toc = time.perf_counter()
                    #timeSpentDoingCxxCall += ( toc - tic)
                    if False:
                        tmpLocalCostAfter, _ = calculator.evaluate_mesh_position( mesh )
                    maximalDeformation = max( maximalDeformation, blockMaximalDeformation )
                    if False:
                        mesh.can_moves = orig_can_moves
                        tmpGlobalCostAfter, _ = calculator.evaluate_mesh_position( mesh )
                        print( f"  block {blockNumber}:" )
                        print( f"    local decrease: {tmpLocalCostBefore-tmpLocalCostAfter} ({tmpLocalCostBefore} - {tmpLocalCostAfter})" )
                        print( f"    global decrease: {tmpGlobalCostBefore-tmpGlobalCostAfter} ({tmpGlobalCostBefore} - {tmpGlobalCostAfter})" )


                globalToc = time.perf_counter()
                print( f"  Total time spent: {globalToc-globalTic:0.4f} sec" )
                #print( f"Total time spent in Cxx call: {timeSpentDoingCxxCall:0.4f} sec ({timeSpentDoingCxxCall/(globalToc-globalTic)*100:0.4f}%)" )

                mesh.can_moves = orig_can_moves
                #print( mesh.can_moves.sum(axis=0) / numberOfNodes * 100 )
                tic = time.perf_counter()
                minLogLikelihoodTimesDeformationPrior, _ = calculator.evaluate_mesh_position( mesh )
                toc = time.perf_counter()
                print( f'  Additional time spent (unnecessarily) computing full cost function: {toc-tic:0.4f} sec' )
                print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                historyOfMaximalDeformation.append(maximalDeformation)
                if maximalDeformation == 0:
                    break

        else:
            print( '=======================================================' )
            print( 'being here' )
            print( '=======================================================' )
          
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
                for blockNumber in range( numberOfBlocks ):
                    tic = time.perf_counter()
                    start = blockNumber * numberOfNodesPerBlock
                    end = min( (blockNumber+1) * numberOfNodesPerBlock, numberOfNodes )
                    
                    mask = np.zeros( numberOfNodes, dtype=bool )
                    mask[ start:end ] = True
                    submesh = mesh.get_submesh( mask );
                    
                    #tmp, _ = calculator.evaluate_mesh_position( mesh )
                    #print( tmp )
                    #tmp, _ = calculator.evaluate_mesh_position( submesh )
                    #print( tmp )
                    
                    #xxx = submesh.points
                    #yyy = mesh.points[ mask, : ]
                    #assert( np.allclose( xxx, yyy ) )
                    #print( xxx.shape )
                    #print( yyy.shape )
                    #print( xxx.max( axis=0 ) )
                    #print( yyy.max( axis=0 ) )
                    #zzz = np.abs( xxx - yyy )
                    #print( zzz.max( axis=0 ) )
                    #print( xxx[ 0,:] )
                    #print( yyy[ 0,:] )
                    #print( xxx[-1,:] )
                    #print( yyy[-1,:] )
                    
                    optimizer = gems.KvlOptimizer(optimizerType, submesh, calculator, optimizationParameters)
                    
                    #submesh.points = mesh.points[ mask, : ]
                    #costxxx, maximalDeformationxxx = optimizer.step_optimizer_samseg()
                    #tmp = mesh.points
                    #tmp[ mask, : ] = submesh.points
                    #mesh.points = tmp
                    #print( costxxx )

                    #submesh.points = mesh.points[ mask, : ]
                    #costxxx, maximalDeformationxxx = optimizer.step_optimizer_samseg()
                    #tmp = mesh.points
                    #tmp[ mask, : ] = submesh.points
                    #mesh.points = tmp
                    #print( costxxx )

                    #submesh.points = mesh.points[ mask, : ]
                    #costxxx, maximalDeformationxxx = optimizer.step_optimizer_samseg()
                    #tmp = mesh.points
                    #tmp[ mask, : ] = submesh.points
                    #mesh.points = tmp
                    #print( costxxx )

                    #import sys
                    #sys.exit()

                    
                    masks.append( mask )
                    submeshes.append( submesh )
                    optimizers.append( optimizer )
                    toc = time.perf_counter()
                    
                    print( f"submesh size {submesh.point_count / numberOfNodesPerBlock * 100:0.4f} %" )
                    print( f"setup time: {toc-tic:0.4f} sec" )
                
                #self.optimizer = optimizers
            else:    
                print( "ProbabilisticAtlas reusing same optimizer!!" )
                optimizers = self.optimizer
                for optimizer in optimizers:
                    optimizer.set_calculator( calculator )
                

            # Run deformation optimization
            historiesOfDeformationCost = [ [] for i in range( numberOfBlocks ) ]
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            
            #
            #numberOfNodes = mesh.point_count
            #mesh.can_moves.sum(axis=0) / numberOfNodes                                                                                          
            #orig_can_moves = mesh.can_moves                                                                                                     
            #numberOfNodesPerBlock = np.ceil( numberOfNodes / numberOfBlocks ).astype( 'int' )
            debug = False
            if debug: debugHistoryOfDeformationCost = []
            currentNodePositions = nodePositionsBeforeDeformation.copy()
            while True:
                #timeSpentDoingCxxCall = 0;
                import time
                globalTic = time.perf_counter()
                maximalDeformation = 0.0
                for blockNumber in range( numberOfBlocks ):
                    if debug:
                        #print( f"blockNumber: {blockNumber}" )
                        mesh.points = currentNodePositions
                        tmpGlobalCostBefore, _ = calculator.evaluate_mesh_position( mesh )

                    submesh = submeshes[ blockNumber ]
                    mask = masks[ blockNumber ]
                    optimizer = optimizers[ blockNumber ]

                    #submesh.points = mesh.points[ mask, : ]
                    submesh.points = currentNodePositions[ mask, : ]
                    

                    if debug:
                        tmpLocalCostBefore, _ = calculator.evaluate_mesh_position( submesh )
                    
                    #tic = time.perf_counter()
                    blockMinLogLikelihoodTimesDeformationPrior, blockMaximalDeformation = optimizer.step_optimizer_samseg()
                    
                    #mesh.points[ mask, : ] = submesh.points
                    #tmp = mesh.points
                    #tmp[ mask, : ] = submesh.points
                    #mesh.points = tmp
                    #xxx = currentNodePositions.copy()
                    currentNodePositions[ mask, : ] = submesh.points
                    #yyy = currentNodePositions.copy()
                    #zzz = np.sum( (xxx-yyy)**2, axis=1 )
                    #print( np.sqrt( zzz.max() ) )
                    #print( blockMaximalDeformation )
                    
                    historiesOfDeformationCost[ blockNumber ].append( blockMinLogLikelihoodTimesDeformationPrior )
                    
                    #toc = time.perf_counter()
                    #timeSpentDoingCxxCall += ( toc - tic)
                    if debug:
                        tmpLocalCostAfter, _ = calculator.evaluate_mesh_position( submesh )
                        #mesh.can_moves = orig_can_moves
                        mesh.points = currentNodePositions
                        tmpGlobalCostAfter, _ = calculator.evaluate_mesh_position( mesh )
                        print( f"  block {blockNumber}:" )
                        print( f"    local decrease: {tmpLocalCostBefore-tmpLocalCostAfter} ({tmpLocalCostBefore} - {tmpLocalCostAfter})" )
                        print( f"    global decrease: {tmpGlobalCostBefore-tmpGlobalCostAfter} ({tmpGlobalCostBefore} - {tmpGlobalCostAfter})" )

                maximalDeformation = max( maximalDeformation, blockMaximalDeformation )
                historyOfMaximalDeformation.append(maximalDeformation)

                globalToc = time.perf_counter()
                print( f"Total time spent: {globalToc-globalTic:0.4f} sec" )
                #print( f"Total time spent in Cxx call: {timeSpentDoingCxxCall:0.4f} sec ({timeSpentDoingCxxCall/(globalToc-globalTic)*100:0.4f}%)" )

                #mesh.can_moves = orig_can_moves
                #print( mesh.can_moves.sum(axis=0) / numberOfNodes * 100 )
                if debug:
                    tic = time.perf_counter()
                    mesh.points = currentNodePositions
                    debugMinLogLikelihoodTimesDeformationPrior, _ = calculator.evaluate_mesh_position( mesh )
                    toc = time.perf_counter()
                    print( f'Additional time spent (unnecessarily) computing full cost function: {toc-tic:0.4f} sec' )
                    print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                    maximalDeformation, debugMinLogLikelihoodTimesDeformationPrior))
                    debugHistoryOfDeformationCost.append(debugMinLogLikelihoodTimesDeformationPrior)
                else:
                    print( f"maximalDeformation={maximalDeformation:.4f}" )                  
                  
                    
                # Check early convergence    
                if maximalDeformation == 0:
                    break
                  
            #      
            mesh.points = currentNodePositions
            minLogLikelihoodTimesDeformationPrior, _ = calculator.evaluate_mesh_position( mesh )
            
            # Reconstruct the historyOfDeformationCost backwards, from incremental improvements in each block
            historyOfDeformationCost = [ minLogLikelihoodTimesDeformationPrior ]
            numberOfIterations = len( historiesOfDeformationCost[0] )
            costOfPreviousIteration = minLogLikelihoodTimesDeformationPrior
            for iterationNumber in reversed( range( 1, numberOfIterations ) ):
                improvementComparedToPreviousIteration = 0.0
                for history in historiesOfDeformationCost:
                    improvementComparedToPreviousIteration += ( history[ iterationNumber-1 ] - history[ iterationNumber ] )
                costOfPreviousIteration += improvementComparedToPreviousIteration
                historyOfDeformationCost.insert( 0, costOfPreviousIteration ) # prepend
        
            if debug:
                assert( np.allclose( np.array( historyOfDeformationCost ), np.array( debugHistoryOfDeformationCost ) ) )
      
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
