import os
import scipy.io
import scipy.ndimage
import numpy as np

import freesurfer as fs
import freesurfer.gems as gems

from .utilities import requireNumpyArray
from .lta import LTA, MRI
from .figures import initVisualizer


def registerAtlas(
        imageFileName,
        meshCollectionFileName,
        templateFileName,
        savePath,
        visualizer=None,
        worldToWorldTransformMatrix=None,
        initLTAFile=None,
        K=None
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
    if visualizer is None: visualizer = initVisualizer(False, False)

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

        # First, we'll hard-code some parameters:
        targetDownsampledVoxelSpacing = 3.0
        # Mesh stiffness K (i.e., measures an average per voxel), so that this needs to be scaled down by the number of voxels that are covered
        if K is None:
            K = 1e-7
        maximalDeformationStopCriterion = 0.005
        lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion  # Doesn't seem to matter very much

        # Initialization
        if initLTAFile is None:
            initialWorldToWorldTransformMatrix = np.identity(4)
        else:
            print('using initial transform: %s' % initLTAFile)
            initLTA = LTA().read(initLTAFile)
            RAS2LPS = np.diag([-1, -1, 1, 1])  # for converting from RAS to LPS
            if initLTA.type == 0:
                # The LTA is vox2vox
                initialWorldToWorldTransformMatrix = RAS2LPS @ initLTA.dstmri.vox2ras0 @ initLTA.xform @ \
                                                     np.linalg.inv(initLTA.srcmri.vox2ras0) @ np.linalg.inv(RAS2LPS)
            else:
                # The LRA is ras2ras
                initialWorldToWorldTransformMatrix = RAS2LPS @ initLTA.xform @ np.linalg.inv(RAS2LPS)

        # Provide an initial (non-identity) affine transform guestimate
        # Rotation around X-axis (direction from left to right ear)
        theta = np.pi / 180 * -10.0
        rotationMatrix = np.identity(4, dtype=np.double)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotationMatrix[1, 1] = cos_theta
        rotationMatrix[1, 2] = -sin_theta
        rotationMatrix[2, 1] = sin_theta
        rotationMatrix[2, 2] = cos_theta
        initialWorldToWorldTransformMatrix = rotationMatrix @ initialWorldToWorldTransformMatrix

        # Isotropic scaling
        scaling = 0.9
        scalingMatrix = np.diag([scaling, scaling, scaling, 1.0])
        initialWorldToWorldTransformMatrix = scalingMatrix @ initialWorldToWorldTransformMatrix

        K = K / (scaling * scaling * scaling)
        multiplied = (initialWorldToWorldTransformMatrix @ templateImageToWorldTransformMatrix)
        initialImageToImageTransformMatrix = np.linalg.solve(imageToWorldTransformMatrix, multiplied)

        # Figure out how much to downsample (depends on voxel size)
        voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)
        downSamplingFactors = np.round(targetDownsampledVoxelSpacing / voxelSpacing)
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
            'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
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

        optimizationSummary = { 'numberOfIterations': len( minLogLikelihoodTimesPriors ), 'cost': minLogLikelihoodTimesPriors[-1] }

    # ------ Save Registration Results ------

    # Save the image-to-image and the world-to-world affine registration matrices
    scipy.io.savemat(os.path.join(savePath, templateFileNameBase + '_coregistrationMatrices.mat'),
                     {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                      'worldToWorldTransformMatrix': worldToWorldTransformMatrix } )

    # Compute and save the talairach.xfm
    lta = computeTalairach(imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix)
    ltaFileName = os.path.join(savePath, 'samseg.talairach.lta')
    print('writing talairach transform to %s' % ltaFileName)
    lta.write(ltaFileName)

    # Save the coregistered template. For historical reasons, we applied the estimated transformation to the template... let's do that now
    desiredTemplateImageToWorldTransformMatrix = np.asfortranarray(imageToWorldTransformMatrix @ imageToImageTransformMatrix)
    transformedTemplateFileName = os.path.join(savePath, templateFileNameBase + '_coregistered' + templateFileNameExtension)
    template.write(transformedTemplateFileName, gems.KvlTransform(desiredTemplateImageToWorldTransformMatrix))

    return worldToWorldTransformMatrix, transformedTemplateFileName, optimizationSummary



def computeTalairach(imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix):
    # Compute the talairach.xfm
    # Load fsaverage orig.mgz -- this is the ultimate target/destination
    fnamedst = os.path.join(fs.fshome(), 'subjects', 'fsaverage', 'mri', 'orig.mgz')
    fsaorig = MRI().read_header(fnamedst)
    # Compute the vox2vox from the template to fsaverage assuming they share world RAS space
    RAS2LPS = np.diag([-1, -1, 1, 1])
    M = np.linalg.inv(RAS2LPS @ fsaorig.vox2ras) @ templateImageToWorldTransformMatrix
    # Compute the input to fsaverage vox2vox by combining the input-template vox2vox and the template-fsaverage vox2vox
    X = M @ np.linalg.inv(imageToImageTransformMatrix)
    # Now write out the LTA. This can be used as the talairach.lta in recon-all
    invol = MRI().read_header(imageFileName)  # have to reread to get header info
    lta = LTA()
    lta.type = 0
    lta.xform = X
    lta.srcfile = imageFileName
    lta.srcmri = invol
    lta.srcmri.vol = []
    lta.dstfile = fnamedst
    lta.dstmri = fsaorig
    lta.dstmri.vol = []
    lta.subject = 'fsaverage'
    return lta
