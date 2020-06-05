import os
import scipy.io
import scipy.ndimage
import numpy as np

import freesurfer as fs
from . import gems, convertLPSTransformToRAS
from .ProbabilisticAtlas import ProbabilisticAtlas

from .utilities import requireNumpyArray
from .figures import initVisualizer
import scipy.ndimage


class Affine:
    def __init__( self, 
                  imageFileName, meshCollectionFileName, templateFileName, 
                  targetDownsampledVoxelSpacing=3.0
                 ):
        self.imageFileName = imageFileName
        self.meshCollectionFileName = meshCollectionFileName
        self.templateFileName = templateFileName
        self.targetDownsampledVoxelSpacing = targetDownsampledVoxelSpacing

    def registerAtlas( self,
            savePath,
            worldToWorldTransformMatrix=None,
            initTransform=None,
            maximalDeformationStopCriterion=0.005,
            Ks=[ 1.0, 0.1 ],
            initializationTryCenterOfGravity=True,
            initializationTryShifts = [ 0, -45, -30, -15, 15, 30, 45 ],
            initializationTryAngles = [ 0 ],
            initializationTryScales = [ 1.0 ],
            initializationTryShiftsSeparately = False,
            visualizer=None ):

        # ------ Set up ------ 
        self.setUp( initTransform, visualizer )

        # ------ Register mesh to image ------
        imageToImageTransformMatrix, worldToWorldTransformMatrix, optimizationSummary = \
            self.registerMeshToImage( worldToWorldTransformMatrix, Ks, maximalDeformationStopCriterion,
                                      initializationTryCenterOfGravity, initializationTryShifts, 
                                      initializationTryAngles, initializationTryScales,
                                      initializationTryShiftsSeparately
                                    )

        # ------ Save results ------
        self.saveResults( savePath, worldToWorldTransformMatrix, imageToImageTransformMatrix )

        return imageToImageTransformMatrix, optimizationSummary
     
     
    #          
    def optimizeTransformation( self, 
                                initialImageToImageTransformMatrix,
                                K, maximalDeformationStopCriterion ):
        
        
          
        # Get the mesh
        initialImageToImageTransform = gems.KvlTransform( requireNumpyArray( initialImageToImageTransformMatrix ) )
        mesh = ProbabilisticAtlas().getMesh( self.meshCollectionFileName, 
                                              transform=initialImageToImageTransform,
                                              K=K )
        originalNodePositions = mesh.points

        
        # Get a registration cost  and stick it in an optimizer
        calculator = gems.KvlCostAndGradientCalculator( 'MutualInformation', [self.image], 'Affine' )
        optimization_parameters = {
            'Verbose': 0,
            'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
            'LineSearchMaximalDeformationIntervalStopCriterion': maximalDeformationStopCriterion,
            'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
        }
        optimizer = gems.KvlOptimizer( 'L-BFGS', mesh, calculator, optimization_parameters)


        # Peform the optimization
        numberOfIterations = 0
        minLogLikelihoodTimesPriors = []
        maximalDeformations = []
        visualizerTitle = f"Affine optimization (K: {K})"
        self.visualizer.start_movie( window_id=visualizerTitle, title=visualizerTitle )
        self.visualizer.show( mesh=mesh, images=self.image.getImageBuffer(), 
                              window_id=visualizerTitle, title=visualizerTitle )
        while True:
            minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_atlas()
            print( "maximalDeformation=%.4f minLogLikelihood=%.4f" % \
                   (maximalDeformation, minLogLikelihoodTimesPrior) )
            minLogLikelihoodTimesPriors.append(minLogLikelihoodTimesPrior)
            maximalDeformations.append(maximalDeformation)

            if maximalDeformation == 0:
                break
            numberOfIterations += 1
            self.visualizer.show( mesh=mesh, images=self.image.getImageBuffer(), 
                                  window_id=visualizerTitle, title=visualizerTitle )

        self.visualizer.show_movie( window_id=visualizerTitle )
        nodePositions = mesh.points
        pointNumbers = [0, 110, 201, 302]
        originalY = np.vstack( ( originalNodePositions[pointNumbers].T, [1, 1, 1, 1] ) )

        Y = np.vstack( ( nodePositions[pointNumbers].T, [1, 1, 1, 1] ) )
        extraImageToImageTransformMatrix = Y @ np.linalg.inv(originalY)
        appliedScaling = np.linalg.det( extraImageToImageTransformMatrix )**(1/3)
        print( f"appliedScaling: {appliedScaling:0.4f}" )
        imageToImageTransformMatrix = extraImageToImageTransformMatrix @ initialImageToImageTransformMatrix 


        optimizationSummary = {'numberOfIterations': len(minLogLikelihoodTimesPriors),
                                'cost': minLogLikelihoodTimesPriors[-1]}
        return imageToImageTransformMatrix, optimizationSummary



    def computeTalairach(self, imageToImageTransformMatrix ):
        # Load fsaverage orig.mgz -- this is the ultimate target/destination
        fnamedst = os.path.join(fs.fshome(), 'subjects', 'fsaverage', 'mri', 'orig.mgz')
        fsaorig = fs.Volume.read(fnamedst)

        # Compute the vox2vox from the template to fsaverage assuming they share world RAS space
        RAS2LPS = np.diag([-1, -1, 1, 1])
        M = np.linalg.inv(RAS2LPS @ fsaorig.affine) @ self.templateImageToWorldTransformMatrix
        # Compute the input to fsaverage vox2vox by combining the input-template vox2vox and the template-fsaverage vox2vox
        vox2vox = M @ np.linalg.inv(imageToImageTransformMatrix)

        # Now write out the LTA. This can be used as the talairach.lta in recon-all
        xform = fs.LinearTransform(vox2vox)
        xform.type = fs.LinearTransform.Type.vox
        xform.source = fs.Volume.read( self.imageFileName ).geometry()
        xform.target = fsaorig.geometry()
        return xform

    #
    def getTransformMatrix( self, angle, scale, shift ):
        #
        # 4x4 matrix mapping [ x' 1 ] to [ y' 1 ] by performing 
        # 
        #    y = S * R * x + t
        # 
        # where R is a 3x3 rotation matrix (rotating around the X-axis -- left-right),
        #       S is a 3x3 isotropic scaling matrix,
        #   and t is a 3x1 translation vector (implementing a shift along the Z-axis -- top-bottom).
        #
        transformMatrix = np.identity( 4, dtype=np.double )
        transformMatrix[ 0, 0 ] = scale
        transformMatrix[ 1, 1 ] = scale * np.cos( angle )
        transformMatrix[ 2, 1 ] = scale * np.sin( angle )
        transformMatrix[ 1, 2 ] = -scale * np.sin( angle )
        transformMatrix[ 2, 2 ] = scale * np.cos( angle )
        transformMatrix[ 2, 3 ] = shift
        
        return transformMatrix
    
    
    #
    def  gridSearch( self, mesh,
                    positionsInTemplateSpace,
                    angles=[ 0.0 ], scales=[ 1.0 ], shifts= [ 0.0 ],
                    baseWorldToWorldTransformMatrices=[ np.eye( 4 ) ],
                    visualizerTitle='Grid search' ):
        #
        # 
        #
        
        # Get a registration cost to evaluate 
        calculator = gems.KvlCostAndGradientCalculator( 'MutualInformation', [ self.image ], 'Affine' )


        #
        bestCost = np.inf
        self.visualizer.start_movie( window_id=visualizerTitle, title=visualizerTitle )
        for angle in angles:
            for scale in scales:
                for baseWorldToWorldTransformMatrix in baseWorldToWorldTransformMatrices:
                    for shift in shifts:
                        # Compute image-to-image transform
                        worldToWorldTransformMatrix = self.getTransformMatrix( angle, scale, shift )
                        imageToImageTransformMatrix = np.linalg.solve( self.imageToWorldTransformMatrix, 
                                                                      worldToWorldTransformMatrix @ \
                                                                      baseWorldToWorldTransformMatrix @ \
                                                                      self.templateImageToWorldTransformMatrix )
      
                        # Apply transform to mesh nodes
                        mesh.points = ProbabilisticAtlas().mapPositionsFromTemplateToSubjectSpace( 
                                        positionsInTemplateSpace, 
                                        gems.KvlTransform( requireNumpyArray( imageToImageTransformMatrix ) ) )
                        self.visualizer.show( images=self.image.getImageBuffer(), mesh=mesh, window_id=visualizerTitle )
                        
                        # Evaluate cost
                        cost, _ = calculator.evaluate_mesh_position( mesh )

                        # If best so far, remember
                        if cost < bestCost:
                            bestCost = cost
                            bestAngle, bestScale, bestShift, bestBaseWorldToWorldTransformMatrix = \
                                angle, scale, shift, baseWorldToWorldTransformMatrix
                            bestImageToImageTransformMatrix = imageToImageTransformMatrix


        # Show the winner
        mesh.points = ProbabilisticAtlas().mapPositionsFromTemplateToSubjectSpace( 
                                    positionsInTemplateSpace, 
                                    gems.KvlTransform( requireNumpyArray( bestImageToImageTransformMatrix ) ) )
        self.visualizer.show( images=self.image.getImageBuffer(), mesh=mesh, window_id=visualizerTitle )
        self.visualizer.show_movie( window_id=visualizerTitle )
        
        #
        return bestAngle, bestScale, bestShift, bestBaseWorldToWorldTransformMatrix, bestImageToImageTransformMatrix
    
    
    #
    def getInitialization( self, tryCenterOfGravity, tryShifts, tryAngles, tryScales, tryShiftsSeparately ):   
        # Get the mesh
        mesh = ProbabilisticAtlas().getMesh( self.meshCollectionFileName, K=0.0 )  # Don't want to use prior on deformation
        positionsInTemplateSpace = mesh.points
        
        baseWorldToWorldTransformMatrices = [ np.eye(4) ]
        if tryCenterOfGravity:
            # Get the centers of gravity of atlas and image, and use that to propose a translation
            templateSize = np.round( np.max( positionsInTemplateSpace, axis=0 ) + 1 ).astype( 'int' )
            priors = mesh.rasterize( templateSize, -1 )
            head = np.sum( priors[:,:,:,1:], axis=3 )
            centerOfGravityTemplate = np.array( scipy.ndimage.measurements.center_of_mass( head ) ) # in image space
            centerOfGravityImage = np.array( 
                scipy.ndimage.measurements.center_of_mass( self.image.getImageBuffer() ) ) # in image space
            centerOfGravityTemplate = self.templateImageToWorldTransformMatrix[ 0:3, 0:3 ] @ centerOfGravityTemplate + \
                                      self.templateImageToWorldTransformMatrix[ 0:3, 3 ] # in world space
            centerOfGravityImage = self.imageToWorldTransformMatrix[ 0:3, 0:3 ] @ centerOfGravityImage + \
                                    self.imageToWorldTransformMatrix[ 0:3, 3 ] # in world space
            centeringWorldToWorldTransformMatrix = np.eye( 4 ); 
            centeringWorldToWorldTransformMatrix[ 0:3, 3 ] = centerOfGravityImage - centerOfGravityTemplate
            baseWorldToWorldTransformMatrices.append( centeringWorldToWorldTransformMatrix )
          
            
        if tryShiftsSeparately:
            # Perform grid search for intialization (done separately for translation and scaling/rotation to limit the
            # number of combinations to be tested)
            _, _, bestShift, bestBaseWorldToWorldTransformMatrix, \
                               _ = self.gridSearch( mesh,
                                                    positionsInTemplateSpace,
                                                    shifts=tryShifts,
                                                    baseWorldToWorldTransformMatrices=baseWorldToWorldTransformMatrices,
                                                    visualizerTitle='Affine grid search shifts' )
              
            _, _, _, _, bestImageToImageTransformMatrix = \
                  self.gridSearch( mesh,
                                   positionsInTemplateSpace,
                                   angles=tryAngles, scales=tryScales, 
                                   shifts = [ bestShift ], 
                                   baseWorldToWorldTransformMatrices=[ bestBaseWorldToWorldTransformMatrix ],
                                   visualizerTitle='Affine grid search anges/scales' )
        else:
            # Basic grid search
            _, _, _, _, bestImageToImageTransformMatrix = \
                  self.gridSearch( mesh,
                                   positionsInTemplateSpace,
                                   angles=tryAngles, scales=tryScales, 
                                   shifts = tryShifts, 
                                   baseWorldToWorldTransformMatrices = baseWorldToWorldTransformMatrices,
                                   visualizerTitle='Affine grid search' )         

        #
        return bestImageToImageTransformMatrix
    

    #
    def setUp( self, initTransform, visualizer ):
        # 
        # Read images etc
        #
      
        # Setup null visualization if necessary
        if visualizer is None:
            visualizer = initVisualizer( False, False )

        # Read in image and template, as well as their coordinates in world (mm) space
        image = gems.KvlImage( self.imageFileName )
        imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
        template = gems.KvlImage( self.templateFileName )
        templateImageToWorldTransformMatrix = template.transform_matrix.as_numpy_array

        # Initialization
        if initTransform is None:
            initialWorldToWorldTransformMatrix = np.identity(4)
        else:
            # Assume the initialization matrix is LPS2LPS
            print('initializing with predifined transform')
            initialWorldToWorldTransformMatrix = initTransform
        originalTemplateImageToWorldTransformMatrix = templateImageToWorldTransformMatrix
        templateImageToWorldTransformMatrix = initialWorldToWorldTransformMatrix @ \
                                                      templateImageToWorldTransformMatrix
        
            
        # Figure out how much to downsample (depends on voxel size)
        voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)
        downSamplingFactors = np.round( self.targetDownsampledVoxelSpacing / voxelSpacing )
        downSamplingFactors[downSamplingFactors < 1] = 1
        upSamplingTranformMatrix = np.diag( np.pad( downSamplingFactors, (0, 1), 
                                                      'constant', constant_values=(0, 1) ) )
        
        # Dowsample image
        imageBuffer = image.getImageBuffer()
        imageBuffer = imageBuffer[ ::int( downSamplingFactors[0] ),
                                  ::int( downSamplingFactors[1] ),
                                  ::int( downSamplingFactors[2] ) ]
        image = gems.KvlImage( requireNumpyArray( imageBuffer ) )
        originalImageToWorldTransformMatrix = imageToWorldTransformMatrix
        imageToWorldTransformMatrix = imageToWorldTransformMatrix @ upSamplingTranformMatrix


        # Save as member variables
        self.visualizer = visualizer
        self.image = image
        self.imageToWorldTransformMatrix = imageToWorldTransformMatrix
        self.templateImageToWorldTransformMatrix = templateImageToWorldTransformMatrix
        self.originalImageToWorldTransformMatrix = originalImageToWorldTransformMatrix
        self.originalTemplateImageToWorldTransformMatrix = originalTemplateImageToWorldTransformMatrix
        self.upSamplingTranformMatrix = upSamplingTranformMatrix
      

    def registerMeshToImage( self, worldToWorldTransformMatrix, Ks, maximalDeformationStopCriterion,
                             tryCenterOfGravity, tryShifts, tryAngles, tryScales, tryShiftsSeparately ):
    
        if worldToWorldTransformMatrix is not None:
            # The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
            # transform (needed for subsequent computations) and be done
            print('world-to-world transform supplied - skipping registration')
            imageToImageTransformMatrix = np.linalg.inv( self.originalImageToWorldTransformMatrix) @ \
                                          worldToWorldTransformMatrix @ \
                                          self.originalTemplateImageToWorldTransformMatrix
            optimizationSummary = None
        else:
            # The solution is not externally given, so we need to compute it.
            print('performing affine atlas registration')
            print('image: %s' % self.imageFileName)
            print('template: %s' % self.templateFileName)

            # Get an initialization of the image-to-image transform (template -> subject)
            import time; tic = time.perf_counter()
            imageToImageTransformMatrix = self.getInitialization( tryCenterOfGravity, tryShifts, 
                                                                  tryAngles, tryScales, tryShiftsSeparately )
            toc = time.perf_counter()
            print( f"Initialization took {toc - tic:0.4f} s" )
            

            # Optimze image-to-image transform (template->subject)
            tic = time.perf_counter()
            for K in Ks:
                imageToImageTransformMatrix, optimizationSummary = \
                    self.optimizeTransformation( imageToImageTransformMatrix,
                                                 K, maximalDeformationStopCriterion )
            toc = time.perf_counter()
            print( f"Optimization took {toc - tic:0.4f} s" )

            
            # Final result: the image-to-image (from template to image) as well as the world-to-world transform that
            # we computed (the latter would be the identity matrix if we didn't move the image at all)
            imageToImageTransformMatrix = self.upSamplingTranformMatrix @ \
                                          imageToImageTransformMatrix # template-to-original-image
            worldToWorldTransformMatrix = self.originalImageToWorldTransformMatrix @ imageToImageTransformMatrix @ \
                                          np.linalg.inv( self.originalTemplateImageToWorldTransformMatrix )
            
        # return result
        return imageToImageTransformMatrix, worldToWorldTransformMatrix, optimizationSummary
      
      
    def  saveResults( self, 
                      savePath, worldToWorldTransformMatrix, imageToImageTransformMatrix ):
                
        # Save the image-to-image and the world-to-world affine registration matrices
        scipy.io.savemat(os.path.join(savePath, 'template_transforms.mat'),
                          {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                          'worldToWorldTransformMatrix': worldToWorldTransformMatrix } )

        # Save the template transform
        inputImage = fs.Volume.read(self.imageFileName)
        templateImage = fs.Volume.read(self.templateFileName)
        xform = fs.LinearTransform(convertLPSTransformToRAS(worldToWorldTransformMatrix))
        xform.type = fs.LinearTransform.Type.ras
        xform.source = templateImage.geometry()
        xform.target = inputImage.geometry()
        ltaFileName = os.path.join(savePath, 'template.lta')
        print('writing template transform to %s' % ltaFileName)
        xform.write(ltaFileName)

        # Compute and save the talairach.xfm
        xform = self.computeTalairach( imageToImageTransformMatrix )
        ltaFileName = os.path.join(savePath, 'samseg.talairach.lta')
        print('writing talairach transform to %s' % ltaFileName)
        xform.write(ltaFileName)

        # Save the coregistered template
        coregistered = fs.Volume(fs.geom.resample(templateImage.data, inputImage.shape[:3], np.linalg.inv(imageToImageTransformMatrix)))
        coregistered.copy_geometry(inputImage)
        coregistered.write(os.path.join(savePath, 'template_coregistered.mgz'))
      
      
