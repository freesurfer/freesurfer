import os
import numpy as np
import scipy.ndimage
import surfa as sf

from freesurfer.samseg import gems
from freesurfer.samseg.warp_mesh import kvlWarpMesh
from freesurfer.samseg.utilities import requireNumpyArray


class ProbabilisticAtlas:
    def __init__(self):
        self.useWarmDeformationOptimizers = False
        self.useBlocksForDeformation = False # Use "grouped coordinate-descent (GCD)" 
                                             # aka "block coordinate-descent" (Fessler)

        if os.environ.get('SAMSEG_EXPERIMENTAL_BLOCK_COORDINATE_DESCENT') is not None: 
            self.useWarmDeformationOptimizers = True
            self.useBlocksForDeformation = True

        self.previousDeformationMesh = None


    def getMesh(self,
                meshCollectionFileName,
                transform=None,
                K=None,
                initialDeformation=None,
                initialDeformationMeshCollectionFileName=None,
                returnInitialDeformationApplied=False,
                competingStructures=None,
                smoothingSigma=0):

        # Get the mesh
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(meshCollectionFileName)
        if K is not None:
            mesh_collection.k = K

        # Do competing structure smoothing if enabled
        if smoothingSigma > 0 and competingStructures:

            print(f'Smoothing competing atlas priors with sigma {smoothingSigma:.2f}')

            # Get initial priors
            size = np.array(mesh_collection.reference_position.max(axis=0) + 1.5, dtype=int)
            priors = mesh_collection.reference_mesh.rasterize(size, -1)

            # Smooth the cortex and WM alphas
            for competingStructureNumbers in competingStructures:
                miniPriors = priors[..., competingStructureNumbers]
                weightsToReassign = np.sum(miniPriors, axis=-1, keepdims=True)
                miniPriors = scipy.ndimage.gaussian_filter(miniPriors.astype(float), sigma=(smoothingSigma, smoothingSigma, smoothingSigma, 0))
                miniPriors /= (np.sum(miniPriors, -1, keepdims=True) + 1e-12)  # Normalize to sum to 1
                miniPriors *= weightsToReassign
                priors[..., competingStructureNumbers] = miniPriors

            # Set the alphas
            alphas = mesh_collection.reference_mesh.fit_alphas(priors, 10)
            mesh_collection.reference_mesh.alphas = alphas

        # Transform
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
            tmp[mask] = data[:, contrastNumber]
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
        
        # Prepare to set up optimizer
        optimizerType = 'L-BFGS'
        optimizationParameters = {
            'Verbose': False,
            'MaximalDeformationStopCriterion': 0.001,  # measured in pixels,
            'LineSearchMaximalDeformationIntervalStopCriterion': 0.001,
            'MaximumNumberOfIterations': 20,
            'BFGS-MaximumMemoryLength': 12
        }
        optimizationParameters.update(userOptimizationParameters)
                
        if not self.useBlocksForDeformation:
            #
            # Vanilla case of a single optimizer
            #
          
            # Get optimizer
            if ( mesh != self.previousDeformationMesh ):
                # Create a new optimizer and plug calculator in it
                optimizer = gems.KvlOptimizer(optimizerType, mesh, calculator, optimizationParameters)
                
                # Remember for a potential next time if warm optimizers are used
                if self.useWarmDeformationOptimizers:
                    self.previousDeformationMesh = mesh
                    self.optimizer = optimizer
                
            else:
                # Reuse a warm optimizer from last time -- just plug the current calculator into it
                print( "ProbabilisticAtlas reusing warm optimizer!!" )
                optimizer = self.optimizer
                optimizer.update_calculator( calculator )
                
            # Run deformation optimization for the specified number of steps (unless it stops moving
            # properly earlier)
            historyOfDeformationCost = []
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            import time
            while True:
                globalTic = time.perf_counter()
                minLogLikelihoodTimesDeformationPrior, maximalDeformation = optimizer.step_optimizer_samseg()
                globalToc = time.perf_counter()
                print("maximalDeformation=%.4f minLogLikelihood=%.4f" % (
                maximalDeformation, minLogLikelihoodTimesDeformationPrior))
                historyOfDeformationCost.append(minLogLikelihoodTimesDeformationPrior)
                print( f"  Total time spent: {globalToc-globalTic:0.4f} sec" )
                historyOfMaximalDeformation.append(maximalDeformation)
                if maximalDeformation == 0:
                    break

        else:
            #
            # Use block coordinate descent in the hope of speeding up convergence.
            # The idea is that solving M optimization problems of size N/M is faster
            # than solving a single big optimization problem of size N. Of course in
            # our case the M problems are not independent, but the coupling is pretty
            # loose (only through a layer of tetrahedra that connect the M point sets)
            #
            numberOfBlocks = 6
            numberOfPassesOverEntireMesh = 3 # Between 1 and optimizationParameters[ 'MaximumNumberOfIterations' ]
          
            # Get optimizers
            if ( mesh != self.previousDeformationMesh ):
                # Create new optimizers and plug calculator in them
                numberOfNodes = mesh.point_count
                numberOfNodesPerBlock = np.ceil( numberOfNodes / numberOfBlocks ).astype( 'int' )
                masks, submeshes, optimizers = [], [], []
                calculators = []
                import time
                useSmartClustering = True
                if useSmartClustering:
                    # Use k-means to group, in an attempt to limit the number of extra points
                    # that need to be added because of boundary tetrahedra connecting 
                    # submeshes with their neighbors
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
                  
                    
                # Remember for a potential next time if warm optimizers are used
                if self.useWarmDeformationOptimizers:
                    self.previousDeformationMesh = mesh
                    self.optimizers = optimizers
                    self.submeshes = submeshes
                    self.masks = masks
 
            else:
                # Reuse warm optimizers from last time -- just plug the current calculator into it
                print( "ProbabilisticAtlas reusing warm optimizers!!" )
                optimizers = self.optimizers
                submeshes = self.submeshes
                masks = self.masks
                for optimizer in optimizers:
                    optimizer.update_calculator( calculator )
                

            # Run deformation optimization for the specified number of steps (unless it stops moving
            # properly earlier)
            historyOfMaximalDeformation = []
            nodePositionsBeforeDeformation = mesh.points
            
            #
            debug = False
            computeHistoryOfDeformationCost = False # Useful for analyzing convergence, but slow
            if computeHistoryOfDeformationCost: historyOfDeformationCost = []
            numberOfInnerLoops = round( optimizationParameters[ 'MaximumNumberOfIterations' ] / 
                                        numberOfPassesOverEntireMesh )
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
                    optimizer.update_mesh( submesh )
                    
                    if debug:
                        tmpLocalCostBefore, _ = calculator.evaluate_mesh_position( submesh )
                    
                  
                    blockPositionsBeforeDeformation = submesh.points
                    for innerLoopNumber in range( numberOfInnerLoops ):
                        optimizer.step_optimizer_samseg()
                    blockPositionsAfterDeformation = submesh.points
                    blockMaximalDeformation = \
                        np.sqrt( np.max( np.sum( ( blockPositionsAfterDeformation - 
                                                   blockPositionsBeforeDeformation ) ** 2, 1 ) ) )
                    
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
                print( f"  Total time spent: {globalToc-globalTic:0.4f} sec" )

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
