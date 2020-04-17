import os
import scipy.io
import scipy.ndimage
import numpy as np

import freesurfer as fs
from . import gems, convertLPSTransformToRAS

from .utilities import requireNumpyArray
from .figures import initVisualizer


class Affine:
    def __init__(self,
                 scaling=0.9,
                 theta=np.pi / 180 * -10.0,
                 K=1e-7,
                 targetDownsampledVoxelSpacing=3.0,
                 maximalDeformationStopCriterion=0.005):
        self.scaling = scaling
        self.theta = theta
        self.targetDownsampledVoxelSpacing = targetDownsampledVoxelSpacing
        self.maximalDeformationStopCriterion = maximalDeformationStopCriterion

    def registerAtlas(self,
            imageFileName,
            meshCollectionFileName,
            templateFileName,
            savePath,
            visualizer=None,
            worldToWorldTransformMatrix=None,
            initTransform=None,
            K=1e-7,
        ):

        # ------ Setup ------

        # Read in image and template, as well as their coordinates in world (mm) space
        image = gems.KvlImage(imageFileName)
        imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
        template = gems.KvlImage(templateFileName)
        templateImageToWorldTransformMatrix = template.transform_matrix.as_numpy_array
        basepath, templateFileNameExtension = os.path.splitext(templateFileName)
        templateFileNameBase = os.path.basename(basepath)

        # Setup null visualization if necessary
        if visualizer is None:
            visualizer = initVisualizer(False, False)

        # ------ Register Image ------

        if worldToWorldTransformMatrix is not None:
            # The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
            # transform (needed for subsequent computations) and be done
            print('world-to-world transform supplied - skipping registration')
            imageToImageTransformMatrix = np.linalg.inv(imageToWorldTransformMatrix) @ worldToWorldTransformMatrix @ templateImageToWorldTransformMatrix
            optimizationSummary = None
        else:
            # The solution is not externally (secretly) given, so we need to compute it.
            print('performing affine atlas registration')
            print('image: %s' % imageFileName)
            print('template: %s' % templateFileName)

            lineSearchMaximalDeformationIntervalStopCriterion = self.maximalDeformationStopCriterion  # Doesn't seem to matter very much

            # Initialization
            if initTransform is None:
                initialWorldToWorldTransformMatrix = np.identity(4)
            else:
                # Assume the initialization matrix is LPS2LPS
                print('initializing with predifined transform')
                initialWorldToWorldTransformMatrix = initTransform

            # Provide an initial (non-identity) affine transform guestestimate
            rotationMatrix, scalingMatrix = self.computeRotationAndScalingMatrixGuessEstimates()

            initialWorldToWorldTransformMatrix = rotationMatrix @ initialWorldToWorldTransformMatrix
            initialWorldToWorldTransformMatrix = scalingMatrix @ initialWorldToWorldTransformMatrix

            K = K / (self.scaling * self.scaling * self.scaling)
            multiplied = (initialWorldToWorldTransformMatrix @ templateImageToWorldTransformMatrix)
            initialImageToImageTransformMatrix = np.linalg.solve(imageToWorldTransformMatrix, multiplied)

            # Figure out how much to downsample (depends on voxel size)
            voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)
            downSamplingFactors = np.round(self.targetDownsampledVoxelSpacing / voxelSpacing)
            downSamplingFactors[downSamplingFactors < 1] = 1

            # Use initial transform to define the reference (rest) position of the mesh (i.e. the one
            # where the log-prior term is zero)
            mesh_collection = gems.KvlMeshCollection()
            mesh_collection.read(meshCollectionFileName)
            mesh_collection.k = K * np.prod(downSamplingFactors)
            mesh_collection.transform(gems.KvlTransform(requireNumpyArray(initialImageToImageTransformMatrix)))
            mesh = mesh_collection.reference_mesh

            # Get image data
            imageBuffer = image.getImageBuffer()
            visualizer.show(images=imageBuffer, window_id='atlas initial', title='Initial Atlas Registration')

            # Downsample
            imageBuffer = imageBuffer[
                          ::int(downSamplingFactors[0]),
                          ::int(downSamplingFactors[1]),
                          ::int(downSamplingFactors[2])]
            image = gems.KvlImage(requireNumpyArray(imageBuffer))
            mesh.scale(1 / downSamplingFactors)
            alphas = mesh.alphas
            gmClassNumber = 3  # Needed for displaying purposes
            numberOfClasses = alphas.shape[1]
            visualizer.show(mesh=mesh, shape=imageBuffer.shape, window_id='atlas mesh', title="Atlas Mesh")

            # Get a registration cost and use it to evaluate some promising starting point proposals
            calculator = gems.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')
            cost, gradient = calculator.evaluate_mesh_position_a(mesh)
            centerOfGravityImage = np.array(scipy.ndimage.measurements.center_of_mass(imageBuffer))
            nodePositions = mesh.points

            # Get the mean node position of the mesh in order to attempt a rough center-of-gravity initialization
            meanNodePosition = np.mean(nodePositions, axis=0)
            baseTranslation = centerOfGravityImage - meanNodePosition

            # Attempt a few different initial alignments. If the initial center-of-gravity placement fails
            # try shifting the translation up and down along the Y axis by 10 voxels
            for verticalDisplacement in (0, -10, 10):
                translation = baseTranslation + np.array((0, verticalDisplacement, 0))
                mesh.points = nodePositions + translation
                trialCost, trialGradient = calculator.evaluate_mesh_position_b(mesh)
                if trialCost >= cost:
                    # Trial alignment was not a success - revert to what we had before
                    mesh.points = nodePositions
                else:
                    # This is better starting position - remember that we applied it
                    initialImageToImageTransformMatrix[0:3, 3] = initialImageToImageTransformMatrix[0:3, 3] + (np.diag(downSamplingFactors) @ translation)
                    break

            originalNodePositions = mesh.points

            # Get an optimizer, and stick the cost function into it
            optimization_parameters = {
                'Verbose': 1.0,
                'MaximalDeformationStopCriterion': self.maximalDeformationStopCriterion,
                'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
                'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
            }
            optimizer = gems.KvlOptimizer( 'L-BFGS', mesh, calculator, optimization_parameters)

            numberOfIterations = 0
            minLogLikelihoodTimesPriors = []
            maximalDeformations = []
            visualizer.start_movie(window_id='atlas iteration', title='Atlas Registration - the movie')
            while True:
                minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_atlas()
                minLogLikelihoodTimesPriors.append(minLogLikelihoodTimesPrior)
                maximalDeformations.append(maximalDeformation)

                if maximalDeformation == 0:
                    break
                numberOfIterations += 1
                visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas iteration', title='Atlas Registration')

            visualizer.show_movie(window_id='atlas iteration')
            nodePositions = mesh.points
            pointNumbers = [0, 110, 201, 302]
            originalY = np.vstack((np.diag(downSamplingFactors) @ originalNodePositions[pointNumbers].T, [1, 1, 1, 1]))

            Y = np.vstack((np.diag(downSamplingFactors) @ nodePositions[pointNumbers].T, [1, 1, 1, 1]))
            extraImageToImageTransformMatrix = Y @ np.linalg.inv(originalY)

            # Final result: the image-to-image (from template to image) as well as the world-to-world transform that
            # we computed (the latter would be the identity matrix if we didn't move the image at all)
            imageToImageTransformMatrix = extraImageToImageTransformMatrix @ initialImageToImageTransformMatrix
            worldToWorldTransformMatrix = imageToWorldTransformMatrix @ imageToImageTransformMatrix @ np.linalg.inv(templateImageToWorldTransformMatrix)

            optimizationSummary = {'numberOfIterations': len(minLogLikelihoodTimesPriors),
                                   'cost': minLogLikelihoodTimesPriors[-1]}

        # ------ Save Registration Results ------

        # Save the image-to-image and the world-to-world affine registration matrices
        scipy.io.savemat(os.path.join(savePath, 'template_transforms.mat'),
                         {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                          'worldToWorldTransformMatrix': worldToWorldTransformMatrix } )

        # Save the template transform
        inputImage = fs.Volume.read(imageFileName)
        templateImage = fs.Volume.read(templateFileName)
        xform = fs.LinearTransform(convertLPSTransformToRAS(worldToWorldTransformMatrix))
        xform.type = fs.LinearTransform.Type.ras
        xform.source = templateImage.geometry()
        xform.target = inputImage.geometry()
        ltaFileName = os.path.join(savePath, 'template.lta')
        print('writing template transform to %s' % ltaFileName)
        xform.write(ltaFileName)

        # Compute and save the talairach.xfm
        xform = self.computeTalairach(imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix)
        ltaFileName = os.path.join(savePath, 'samseg.talairach.lta')
        print('writing talairach transform to %s' % ltaFileName)
        xform.write(ltaFileName)

        # Save the coregistered template
        coregistered = fs.Volume(fs.geom.resample(templateImage.data, inputImage.shape[:3], np.linalg.inv(imageToImageTransformMatrix)))
        coregistered.copy_geometry(inputImage)
        coregistered.write(os.path.join(savePath, 'template_coregistered.mgz'))

        # Historically, the template was resaved with a transformed header and the image-to-image transform
        # was later extracted by comparing the input image and coregistered template transforms, but shear
        # cannot be saved through an ITK image, so a better method is to pass the image-to-image transform matrix
        # directly to samseg. If the SAMSEG_LEGACY_REGISTRATION env var is set, this old method is enabled
        if os.environ.get('SAMSEG_LEGACY_REGISTRATION') is not None:
            print('INFO: using legacy (broken) registration option')
            desiredTemplateImageToWorldTransformMatrix = np.asfortranarray(imageToWorldTransformMatrix @ imageToImageTransformMatrix)
            transformedTemplateFileName = os.path.join(savePath, 'template_coregistered_legacy.nii')
            template.write(transformedTemplateFileName, gems.KvlTransform(desiredTemplateImageToWorldTransformMatrix))

        return imageToImageTransformMatrix, optimizationSummary

    def computeRotationAndScalingMatrixGuessEstimates(self):
        # Rotation around X-axis (direction from left to right ear)
        rotationMatrix = np.identity(4, dtype=np.double)
        cos_theta = np.cos(self.theta)
        sin_theta = np.sin(self.theta)
        rotationMatrix[1, 1] = cos_theta
        rotationMatrix[1, 2] = -sin_theta
        rotationMatrix[2, 1] = sin_theta
        rotationMatrix[2, 2] = cos_theta

        # Isotropic scaling
        scalingMatrix = np.diag([self.scaling, self.scaling, self.scaling, 1.0])

        return rotationMatrix, scalingMatrix

    def computeTalairach(self, imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix):
        # Load fsaverage orig.mgz -- this is the ultimate target/destination
        fnamedst = os.path.join(fs.fshome(), 'subjects', 'fsaverage', 'mri', 'orig.mgz')
        fsaorig = fs.Volume.read(fnamedst)

        # Compute the vox2vox from the template to fsaverage assuming they share world RAS space
        RAS2LPS = np.diag([-1, -1, 1, 1])
        M = np.linalg.inv(RAS2LPS @ fsaorig.affine) @ templateImageToWorldTransformMatrix
        # Compute the input to fsaverage vox2vox by combining the input-template vox2vox and the template-fsaverage vox2vox
        vox2vox = M @ np.linalg.inv(imageToImageTransformMatrix)

        # Now write out the LTA. This can be used as the talairach.lta in recon-all
        xform = fs.LinearTransform(vox2vox)
        xform.type = fs.LinearTransform.Type.vox
        xform.source = fs.Volume.read(imageFileName).geometry()
        xform.target = fsaorig.geometry()
        return xform
