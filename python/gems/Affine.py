import os
import scipy.io
import scipy.ndimage
import numpy as np
import surfa as sf

import gems
from gems import convertLPSTransformToRAS
from gems.ProbabilisticAtlas import ProbabilisticAtlas
from gems.utilities import requireNumpyArray
from gems.figures import initVisualizer


class initializationOptions:
    def __init__( self,
                  pitchAngles=np.array( [ 0 ] ) / 180.0 * np.pi, # in radians
                  scales=[ 1.0 ],  
                  horizontalTableShifts=[ 40.0, 20.0, 0.0, -20.0, -40.0 ], # in mm, in template space
                  verticalTableShifts=[ 0.0 ], # in mm, in template space
                  tryCenterOfGravity=True,
                  searchForTableShiftsSeparately=False,
                  pitchCenter=[ 0.0, 0.0, 0.0 ], # in mm, in template space - anterior commissure
                  scalingCenter=[ 0.0, 120.0, 0.0 ], # in mm, in template space - back of the head
                  initialPitchAngle=-10.0/180.0*np.pi, # in radians
                  initialScale=0.9, 
                  initialTableShift=[ 0.0, 0.0, 0.0 ] # in mm, in template space
                ):
        self.pitchAngles = pitchAngles
        self.scales = scales
        self.horizontalTableShifts = horizontalTableShifts
        self.verticalTableShifts = verticalTableShifts
        self.tryCenterOfGravity = tryCenterOfGravity
        self.searchForTableShiftsSeparately = searchForTableShiftsSeparately
        self.pitchCenter = pitchCenter
        self.scalingCenter = scalingCenter
        self.initialPitchAngle = initialPitchAngle
        self.initialScale = initialScale
        self.initialTableShift = initialTableShift




class Affine:
    def __init__( self, 
                  imageFileName, meshCollectionFileName, templateFileName ):
        self.imageFileName = imageFileName
        self.meshCollectionFileName = meshCollectionFileName
        self.templateFileName = templateFileName

    def registerAtlas( self,
            savePath,
            worldToWorldTransformMatrix=None,
            initTransform=None,
            Ks=[ 20.0, 10.0, 5.0 ],
            initializationOptions = initializationOptions(),
            targetDownsampledVoxelSpacing=3.0,
            maximalDeformationStopCriterion=0.005,
            visualizer=None ):

        # ------ Set up ------ 
        self.setUp( initTransform, targetDownsampledVoxelSpacing, visualizer )

        # ------ Register mesh to image ------
        imageToImageTransformMatrix, worldToWorldTransformMatrix, optimizationSummary = \
            self.registerMeshToImage( worldToWorldTransformMatrix=worldToWorldTransformMatrix, 
                                      Ks=Ks, 
                                      maximalDeformationStopCriterion=maximalDeformationStopCriterion,
                                      initializationOptions=initializationOptions )

        # ------ Save results ------
        self.saveResults( savePath, worldToWorldTransformMatrix, imageToImageTransformMatrix )

        return imageToImageTransformMatrix, optimizationSummary
     
     
    #          
    def optimizeTransformation( self, 
                                initialImageToImageTransformMatrix,
                                K, maximalDeformationStopCriterion ):
        
        # In our models the mesh stiffness (scaling the log-prior) is relative to the log-likelihood,
        # which scales with the number of voxels being modeled. But for MI the log-likelihood
        # is normalized, so we need to divide the mesh stiffness by the (expected) number of voxels
        # covered in order to compensate.
        Keffective = K / self.expectedNumberOfVoxelsCovered
          
        # Get the mesh
        initialImageToImageTransform = gems.KvlTransform( requireNumpyArray( initialImageToImageTransformMatrix ) )
        mesh = ProbabilisticAtlas().getMesh( self.meshCollectionFileName, 
                                              transform=initialImageToImageTransform,
                                              K=Keffective )
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
        fnamedst = os.path.join(os.environ.get('FREESURFER_HOME'), 'subjects', 'fsaverage', 'mri', 'orig.mgz')
        fsaorig = sf.load_volume(fnamedst)

        # Compute the vox2vox from the template to fsaverage assuming they share world RAS space
        RAS2LPS = np.diag([-1, -1, 1, 1])
        M = np.linalg.inv(RAS2LPS @ fsaorig.geom.vox2world) @ self.templateImageToWorldTransformMatrix
        # Compute the input to fsaverage vox2vox by combining the input-template vox2vox and the template-fsaverage vox2vox
        vox2vox = M @ np.linalg.inv(imageToImageTransformMatrix)

        # Now write out the LTA. This can be used as the talairach.lta in recon-all
        return sf.Affine(vox2vox, space='vox', source=sf.load_volume(self.imageFileName), target=fsaorig)

    def getTransformMatrix( self, 
                            pitchAngle=0.0, scale=1.0, horizontalTableShift=0.0, verticalTableShift=0.0,
                            initialTableShift=None ):
        #
        # 4x4 matrix simulating translation (horizontal and vertical table shift of the MRI scanner), 
        # followed by isotropic scaling (around scalingCenter), and finally rotation (pitch around pitchCenter)
        # of the atlas template. Everything is measured in the template world coordinates (i.e., measured in mm).
        #       
        if initialTableShift is None:
            initialTableShift = self.initialTableShift
        
        
        # Translation
        translationMatrix = np.identity( 4, dtype=np.double )
        translationMatrix[ 0:3, 3 ] = initialTableShift
        translationMatrix[ 1, 3 ] += verticalTableShift
        translationMatrix[ 2, 3 ] += horizontalTableShift
        
        # Scaling
        scalingMatrix = np.identity( 4, dtype=np.double )
        scalingMatrix[ 0, 0 ] = scale * self.initialScale
        scalingMatrix[ 1, 1 ] = scale * self.initialScale
        scalingMatrix[ 2, 2 ] = scale * self.initialScale
        scalingCenter = np.array( [ self.scalingCenter[0], self.scalingCenter[1], self.scalingCenter[2], 
                                    1 ] ).reshape( -1, 1 )
        scalingCenter = translationMatrix @ scalingCenter
        centeringMatrix = np.eye( 4 )
        centeringMatrix[ 0:3, 3 ] = -scalingCenter[ 0:3, 0 ]
        scalingMatrix = np.linalg.solve( centeringMatrix, scalingMatrix @ centeringMatrix )
            
        # Rotation
        rotationMatrix = np.identity( 4, dtype=np.double )
        rotationMatrix[ 1, 1 ] = np.cos( pitchAngle + self.initialPitchAngle )
        rotationMatrix[ 2, 1 ] = np.sin( pitchAngle + self.initialPitchAngle )
        rotationMatrix[ 1, 2 ] = -np.sin( pitchAngle + self.initialPitchAngle )
        rotationMatrix[ 2, 2 ] = np.cos( pitchAngle + self.initialPitchAngle )
        rotationCenter = np.array( [ self.pitchCenter[0], self.pitchCenter[1], self.pitchCenter[2], 
                                     1 ] ).reshape( -1, 1 )
        rotationCenter = scalingMatrix @ translationMatrix @ rotationCenter
        centeringMatrix = np.eye( 4 )
        centeringMatrix[ 0:3, 3 ] = -rotationCenter[ 0:3, 0 ]
        rotationMatrix = np.linalg.solve( centeringMatrix, rotationMatrix @ centeringMatrix )
        

        transformMatrix = rotationMatrix @ scalingMatrix @ translationMatrix
        
        return transformMatrix
    
    
    #
    def gridSearch( self, mesh,
                    positionsInTemplateSpace,
                    pitchAngles=[ 0.0 ],
                    scales= [ 1.0 ],
                    horizontalTableShifts = [ 0.0 ],
                    verticalTableShifts = [ 0.0 ],
                    initialTableShifts=[ 0.0, 0.0, 0.0 ],
                    visualizerTitle='Grid search' ):
        #
        # 
        #
        
        # Get a registration cost to evaluate 
        calculator = gems.KvlCostAndGradientCalculator( 'MutualInformation', [ self.image ], 'Affine' )


        #
        bestCost = np.inf
        self.visualizer.start_movie( window_id=visualizerTitle, title=visualizerTitle )
        for initialTableShift in initialTableShifts:
            for scale in scales:
                for horizontalTableShift in horizontalTableShifts:
                    for verticalTableShift in verticalTableShifts:
                        for pitchAngle in pitchAngles:
                            # Compute image-to-image transform
                            worldToWorldTransformMatrix = \
                                self.getTransformMatrix( pitchAngle=pitchAngle, 
                                                         scale=scale, 
                                                         horizontalTableShift=horizontalTableShift,
                                                         verticalTableShift=verticalTableShift,
                                                         initialTableShift=initialTableShift )
                            imageToImageTransformMatrix = np.linalg.solve( self.imageToWorldTransformMatrix,
                                                                           worldToWorldTransformMatrix @ \
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
                                bestPitchAngle, bestScale, bestHorizontalTableShift, \
                                  bestVerticalTableShift, bestInitialTableShift = \
                                    pitchAngle, scale, horizontalTableShift, \
                                      verticalTableShift, initialTableShift
                                bestImageToImageTransformMatrix = imageToImageTransformMatrix


        # Show the winner
        mesh.points = ProbabilisticAtlas().mapPositionsFromTemplateToSubjectSpace( 
                                    positionsInTemplateSpace, 
                                    gems.KvlTransform( requireNumpyArray( bestImageToImageTransformMatrix ) ) )
        self.visualizer.show( images=self.image.getImageBuffer(), mesh=mesh, window_id=visualizerTitle )
        self.visualizer.show_movie( window_id=visualizerTitle )
        
        #
        return bestPitchAngle, bestScale, bestHorizontalTableShift, bestVerticalTableShift, \
                 bestInitialTableShift, bestImageToImageTransformMatrix
    
    
    #
    def getInitialization( self, initializationOptions ):
        # Remember some things that are fixed
        self.pitchCenter = initializationOptions.pitchCenter
        self.scalingCenter = initializationOptions.scalingCenter
        self.initialTableShift = initializationOptions.initialTableShift
        self.initialPitchAngle = initializationOptions.initialPitchAngle
        self.initialScale = initializationOptions.initialScale
      
      
        # Get the mesh. Although we're not penalizing prior deformation here (K=0), we still need to first put 
        # the reference mesh in subject space to make sure no tetrahedra get negative volume (which would trigger
        # a sentinel "infinity" cost even when K=0)
        imageToImageTransformMatrix = np.linalg.inv( self.imageToWorldTransformMatrix) @ \
                                      self.templateImageToWorldTransformMatrix
        imageToImageTransform = gems.KvlTransform( requireNumpyArray( imageToImageTransformMatrix ) )
        mesh = ProbabilisticAtlas().getMesh( self.meshCollectionFileName, 
                                             transform=imageToImageTransform,
                                             K=0.0 )
        positionsInTemplateSpace = ProbabilisticAtlas().mapPositionsFromSubjectToTemplateSpace( mesh.points, 
                                                                          imageToImageTransform )

        
        #
        bestInitialTableShift = initializationOptions.initialTableShift
        if initializationOptions.tryCenterOfGravity:
            # Get the centers of gravity of atlas and image, and use that to propose a translation
            mesh.points = positionsInTemplateSpace
            templateSize = np.round( np.max( mesh.points, axis=0 ) + 1 ).astype( 'int' )
            priors = mesh.rasterize( templateSize, -1 )
            head = np.sum( priors[:,:,:,1:], axis=3 )
            centerOfGravityTemplate = np.array( scipy.ndimage.measurements.center_of_mass( head ) ) # in image space
            centerOfGravityImage = np.array( 
                scipy.ndimage.measurements.center_of_mass( self.image.getImageBuffer() ) ) # in image space
            tmp = self.getTransformMatrix()
            tmp = tmp @ self.templateImageToWorldTransformMatrix
            centerOfGravityTemplate = tmp[ 0:3, 0:3 ] @ centerOfGravityTemplate + \
                                      tmp[ 0:3, 3 ] # in world space
            centerOfGravityImage = self.imageToWorldTransformMatrix[ 0:3, 0:3 ] @ centerOfGravityImage + \
                                   self.imageToWorldTransformMatrix[ 0:3, 3 ] # in world space
            centeringTableShift = centerOfGravityImage - centerOfGravityTemplate
                
            
            # Try which one is best
            initialTableShifts = [ initializationOptions.initialTableShift, centeringTableShift ]
            _, _, _, _, bestInitialTableShift, _ = \
                                    self.gridSearch( mesh,
                                                     positionsInTemplateSpace,
                                                     initialTableShifts=initialTableShifts,
                                                     visualizerTitle='Affine grid search center of gravity' )
            
            
        if initializationOptions.searchForTableShiftsSeparately:
            # Perform grid search for intialization (done separately for translation and scaling/rotation to limit the
            # number of combinations to be tested)
            _, _, bestHorizontalTableShift, bestVerticalTableShift, _, _ = \
                  self.gridSearch( mesh,
                                   positionsInTemplateSpace,
                                   horizontalTableShifts=initializationOptions.horizontalTableShifts,
                                   verticalTableShifts=initializationOptions.verticalTableShifts,
                                   initialTableShifts=[ bestInitialTableShift ],
                                   visualizerTitle='Affine grid search shifts' )

            _, _, _, _, _, bestImageToImageTransformMatrix = \
                  self.gridSearch( mesh,
                                   positionsInTemplateSpace,
                                   pitchAngles=initializationOptions.pitchAngles,
                                   scales=initializationOptions.scales,
                                   horizontalTableShifts = [ bestHorizontalTableShift ],
                                   verticalTableShifts = [ bestVerticalTableShift ],
                                   initialTableShifts=[ bestInitialTableShift ],
                                   visualizerTitle='Affine grid search anges/scales' )
        else:
            # Basic grid search
            _, _, _, _, _, bestImageToImageTransformMatrix = \
                  self.gridSearch( mesh,
                                   positionsInTemplateSpace,
                                   pitchAngles=initializationOptions.pitchAngles,
                                   scales=initializationOptions.scales,
                                   horizontalTableShifts=initializationOptions.horizontalTableShifts,
                                   verticalTableShifts=initializationOptions.verticalTableShifts,
                                   initialTableShifts=[ bestInitialTableShift ],
                                   visualizerTitle='Affine grid search' )         

        #
        return bestImageToImageTransformMatrix
    

    #
    def setUp( self, initTransform, targetDownsampledVoxelSpacing, visualizer ):
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

        # Initialization -- since we're attaching physical meaning to initializations in template space,
        # let's not touch the template but rather move the subject image according to the initialization
        if initTransform is None:
            initialWorldToWorldTransformMatrix = np.identity(4)
        else:
            # Assume the initialization matrix is LPS2LPS
            print('initializing with predifined transform')
            initialWorldToWorldTransformMatrix = initTransform
        originalImageToWorldTransformMatrix = imageToWorldTransformMatrix
        imageToWorldTransformMatrix = np.linalg.solve( initialWorldToWorldTransformMatrix, 
                                                       imageToWorldTransformMatrix )                                           
        
            
        # Figure out how much to downsample (depends on voxel size)
        voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)
        downSamplingFactors = np.round( targetDownsampledVoxelSpacing / voxelSpacing )
        downSamplingFactors[downSamplingFactors < 1] = 1
        upSamplingTranformMatrix = np.diag( np.pad( downSamplingFactors, (0, 1), 
                                                      'constant', constant_values=(0, 1) ) )
        
        # Dowsample image
        imageBuffer = image.getImageBuffer()
        imageBuffer = imageBuffer[ ::int( downSamplingFactors[0] ),
                                  ::int( downSamplingFactors[1] ),
                                  ::int( downSamplingFactors[2] ) ]
        image = gems.KvlImage( requireNumpyArray( imageBuffer ) )
        imageToWorldTransformMatrix = imageToWorldTransformMatrix @ upSamplingTranformMatrix


        # Save as member variables
        self.visualizer = visualizer
        self.image = image
        self.imageToWorldTransformMatrix = imageToWorldTransformMatrix
        self.templateImageToWorldTransformMatrix = templateImageToWorldTransformMatrix
        self.originalImageToWorldTransformMatrix = originalImageToWorldTransformMatrix
        self.upSamplingTranformMatrix = upSamplingTranformMatrix
        self.expectedNumberOfVoxelsCovered = np.prod( template.getImageBuffer().shape )
         
         

    def registerMeshToImage( self, worldToWorldTransformMatrix, Ks, maximalDeformationStopCriterion,
                             initializationOptions ):
    
        if worldToWorldTransformMatrix is not None:
            # The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
            # transform (needed for subsequent computations) and be done
            print('world-to-world transform supplied - skipping registration')
            imageToImageTransformMatrix = np.linalg.inv( self.originalImageToWorldTransformMatrix) @ \
                                          worldToWorldTransformMatrix @ \
                                          self.templateImageToWorldTransformMatrix
            optimizationSummary = None
        else:
            # The solution is not externally given, so we need to compute it.
            print('performing affine atlas registration')
            print('image: %s' % self.imageFileName)
            print('template: %s' % self.templateFileName)

            # Get an initialization of the image-to-image transform (template -> subject)
            imageToImageTransformMatrix = self.getInitialization( initializationOptions )


            # Optimze image-to-image transform (template->subject)
            for K in Ks:
                imageToImageTransformMatrix, optimizationSummary = \
                    self.optimizeTransformation( imageToImageTransformMatrix,
                                                 K, maximalDeformationStopCriterion )

            
            # Final result: the image-to-image (from template to image) as well as the world-to-world transform that
            # we computed (the latter would be the identity matrix if we didn't move the image at all)
            imageToImageTransformMatrix = self.upSamplingTranformMatrix @ \
                                          imageToImageTransformMatrix # template-to-original-image
            worldToWorldTransformMatrix = self.originalImageToWorldTransformMatrix @ imageToImageTransformMatrix @ \
                                          np.linalg.inv( self.templateImageToWorldTransformMatrix )
            
        # return result
        return imageToImageTransformMatrix, worldToWorldTransformMatrix, optimizationSummary
      
      
    def saveResults(self, savePath, worldToWorldTransformMatrix, imageToImageTransformMatrix):
                
        # Save the image-to-image and the world-to-world affine registration matrices
        scipy.io.savemat(os.path.join(savePath, 'template_transforms.mat'),
                          {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                          'worldToWorldTransformMatrix': worldToWorldTransformMatrix } )

        # Save the template transform
        inputImage = sf.load_volume(self.imageFileName)
        templateImage = sf.load_volume(self.templateFileName)
        xform = sf.Affine(
            convertLPSTransformToRAS(worldToWorldTransformMatrix),
            space='world',
            source=templateImage,
            target=inputImage,
        )
        ltaFileName = os.path.join(savePath, 'template.lta')
        print('writing template transform to %s' % ltaFileName)
        xform.save(ltaFileName)

        # Compute and save the talairach.xfm
        xform = self.computeTalairach( imageToImageTransformMatrix )
        ltaFileName = os.path.join(savePath, 'samseg.talairach.lta')
        print('writing talairach transform to %s' % ltaFileName)
        xform.save(ltaFileName)

        # Save the coregistered template
        xform = sf.Affine(imageToImageTransformMatrix, space='voxel', source=templateImage, target=inputImage)
        coregistered = templateImage.transform(xform)
        coregistered.save(os.path.join(savePath, 'template_coregistered.mgz'))
