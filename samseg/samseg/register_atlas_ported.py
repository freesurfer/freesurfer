import logging
from gems2python import GEMS2Python

import numpy as np
import scipy.ndimage
import scipy.io
import os

from samseg.lta import LTA, MRI
from samseg.show_figures import DoNotShowFigures

logger = logging.getLogger(__name__)


def require_np_array(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])


ASSERTIONS_ON = False

SKIP_ATLAS_SHOW_FIGURES = False

def assert_close(golden, trial, **kwargs):
    if ASSERTIONS_ON:
        np.testing.assert_allclose(golden, trial, **kwargs)


def samseg_registerAtlas(imageFileName,
                         meshCollectionFileName,
                         templateFileName,
                         savePath,
                         visualizer=None,
                         worldToWorldTransformMatrix=None,
                         InitLTAFile=None,
                         checkpoint_manager=None):
    if SKIP_ATLAS_SHOW_FIGURES or visualizer is None:
        visualizer = DoNotShowFigures()

    # For converting from RAS to LPS. Itk/GEMS/SAMSEG uses LPS internally
    RAS2LPS = np.diag([-1, -1, 1, 1])

    # Print out the input
    print('entering registerAtlas')
    print(imageFileName)
    print(meshCollectionFileName)
    print(templateFileName)
    print(savePath)
    print(visualizer)
    print(worldToWorldTransformMatrix)
    if InitLTAFile:
        print(InitLTAFile)
    else:
        print('No init.lta file')

    # Read in image and template, as well as their coordinates in world (mm) space
    image = GEMS2Python.KvlImage(imageFileName)
    imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
    template = GEMS2Python.KvlImage(templateFileName)

    templateImageToWorldTransformMatrix = template.transform_matrix.as_numpy_array
    basepath, templateFileNameExtension = os.path.splitext(templateFileName)
    _, templateFileNameBase = os.path.split(basepath)
    if worldToWorldTransformMatrix is None:
        # The solution is not externally (secretly) given, so we need to compute it
        #
        #
        # Some hard-coded parameter settings first
        targetDownsampledVoxelSpacing = 3.0
        K = 1e-7 # Mesh stiffness
        #          (i.e., measures an average *per voxel*), so that this needs to be scaled down by the
        #          number of voxels that are covered
        maximalDeformationStopCriterion = 0.005
        lineSearchMaximalDeformationIntervalStopCriterion = maximalDeformationStopCriterion # Doesn't seem to matter very much
        #
        #
        # Initialization
        if InitLTAFile is None:
            initialWorldToWorldTransformMatrix = np.identity(4)
        else:
            InitLTA = LTA().read(InitLTAFile)
            if InitLTA.type == 0:
                # The LTA is vox2vox
                initialWorldToWorldTransformMatrix = \
                    RAS2LPS @ InitLTA.dstmri.vox2ras0 @ \
                    InitLTA.xform @ \
                    np.linalg.inv(InitLTA.srcmri.vox2ras0) @ np.linalg.inv(RAS2LPS)
            else:
                # The LRA is ras2ras
                initialWorldToWorldTransformMatrix = \
                    RAS2LPS @ \
                    InitLTA.xform @ \
                    np.linalg.inv(RAS2LPS)
        #     % Provide an initial (non-identity) affine transform guestimate
        #     % Rotation around X-axis (direction from left to right ear)
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
        # Use initial transform to define the reference (rest) position of the mesh (i.e., the one
        # where the log-prior term is zero)
        mesh_collection = GEMS2Python.KvlMeshCollection()
        mesh_collection.read(meshCollectionFileName)
        mesh_collection.k = K * np.prod(downSamplingFactors)
        mesh_collection.transform(GEMS2Python.KvlTransform(require_np_array(initialImageToImageTransformMatrix)))
        mesh = mesh_collection.reference_mesh
        #   % Get image data
        imageBuffer = image.getImageBuffer()
        visualizer.show(images=imageBuffer, window_id='atlas initial', title='Initial Atlas Registration')
        #   % Downsample
        imageBuffer = imageBuffer[
                      ::int(downSamplingFactors[0]),
                      ::int(downSamplingFactors[1]),
                      ::int(downSamplingFactors[2])]
        image = GEMS2Python.KvlImage(require_np_array(imageBuffer))
        mesh.scale(1 / downSamplingFactors)
        alphas = mesh.alphas
        gmClassNumber = 3  # Needed for displaying purposes
        numberOfClasses = alphas.shape[1]
        visualizer.show(mesh=mesh, shape=imageBuffer.shape, window_id='atlas mesh', title="Atlas Mesh")
        # Get a registration cost and use it to evaluate some promising starting point proposals
        calculator = GEMS2Python.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')
        cost, gradient = calculator.evaluate_mesh_position_a(mesh)
        centerOfGravityImage = np.array(scipy.ndimage.measurements.center_of_mass(imageBuffer))
        priors = mesh.rasterize_atlas(imageBuffer.shape)
        visualizer.show(probabilities=priors, window_id='atlas probabilities', title='Atlas Probabilities')
        tmp = np.sum(priors[:, :, :, 1:], axis=3)
        centerOfGravityAtlas = np.array(scipy.ndimage.measurements.center_of_mass(tmp))
        initialTranslation = centerOfGravityImage - centerOfGravityAtlas
        nodePositions = mesh.points
        trialNodePositions = nodePositions + initialTranslation
        mesh.points = trialNodePositions
        trialCost, trialGradient = calculator.evaluate_mesh_position_b(mesh)
        if trialCost >= cost:
            # Center of gravity was not a success; revert to what we had before
            mesh.points = nodePositions
        else:
            # This is better starting position; remember that we applied it
            initialImageToImageTransformMatrix[0:3, 3] = initialImageToImageTransformMatrix[0:3, 3] + (
                    np.diag(downSamplingFactors) @ initialTranslation)
        originalNodePositions = mesh.points

        # Get an optimizer, and stick the cost function into it
        optimization_parameters = {
            'Verbose': 1.0,
            'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
            'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
            'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
        }
        optimizer = GEMS2Python.KvlOptimizer(
            'L-BFGS',
            mesh,
            calculator,
            optimization_parameters
        )

        numberOfIterations = 0
        minLogLikelihoodTimesPriors = []
        maximalDeformations = []
        costs = []
        gradients = []
        visualizer.start_movie(window_id='atlas iteration', title='Atlas Registration - the movie')
        while True:
            cost, gradient = calculator.evaluate_mesh_position_c(mesh)
            logger.info("cost = %f", cost)
            costs.append(cost)
            gradients.append(gradient)

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
        worldToWorldTransformMatrix = imageToWorldTransformMatrix @ imageToImageTransformMatrix @ np.linalg.inv(
            templateImageToWorldTransformMatrix)

    else:
        # The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
        # transform (needed for subsequent computations) and be done
        imageToImageTransformMatrix = np.linalg.inv(
            imageToWorldTransformMatrix) * worldToWorldTransformMatrix @ templateImageToWorldTransformMatrix

    transformedTemplateFileName = save_results(
        costs,
        imageFileName,
        imageToImageTransformMatrix,
        imageToWorldTransformMatrix,
        savePath,
        template,
        templateFileNameBase,
        templateFileNameExtension,
        templateImageToWorldTransformMatrix,
        worldToWorldTransformMatrix,
    )
    return worldToWorldTransformMatrix, transformedTemplateFileName


def save_results(
        costs,
        imageFileName,
        imageToImageTransformMatrix,
        imageToWorldTransformMatrix,
        savePath,
        template,
        templateFileNameBase,
        templateFileNameExtension,
        templateImageToWorldTransformMatrix,
        worldToWorldTransformMatrix,
):
    save_coregistration_matrices(
        costs,
        imageToImageTransformMatrix,
        savePath,
        templateFileNameBase,
        worldToWorldTransformMatrix,
    )
    save_talairch(
        compute_talairach(imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix),
        savePath)
    transformedTemplateFileName = save_coregistered_template(
        imageToImageTransformMatrix,
        imageToWorldTransformMatrix,
        savePath,
        template,
        templateFileNameBase,
        templateFileNameExtension
    )
    return transformedTemplateFileName


def save_coregistration_matrices(costs, imageToImageTransformMatrix, savePath, templateFileNameBase,
                                 worldToWorldTransformMatrix):
    # Save the image-to-image and the world-to-world affine registration matrices
    scipy.io.savemat(os.path.join(savePath, templateFileNameBase + '_coregistrationMatrices.mat'),
                     {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                      'worldToWorldTransformMatrix': worldToWorldTransformMatrix,
                      'costs': costs,
                      }
                     )


def compute_talairach(imageFileName, imageToImageTransformMatrix, templateImageToWorldTransformMatrix, fshome=None):
    # Compute the talairach.xfm
    # Load fsaverage orig.mgz -- this is the ultimate target/destination
    if fshome is None:
        fshome = os.getenv('FREESURFER_HOME')
    fnamedst = os.path.join(fshome, 'subjects', 'fsaverage', 'mri', 'orig.mgz')
    fsaorig = MRI().read_header(fnamedst)
    # Compute the vox2vox from the template to fsaverage assuming they
    # share world RAS space
    RAS2LPS = np.diag([-1, -1, 1, 1])
    M = np.linalg.inv(RAS2LPS @ fsaorig.vox2ras) @ templateImageToWorldTransformMatrix

    # Compute the input to fsaverage vox2vox by combining the
    # input-template vox2vox and the template-fsaverage vox2vox
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


def save_talairch(lta, savePath):
    ltaFileName = os.path.join(savePath, 'samseg.talairach.lta')
    lta.write(ltaFileName)
    print('Done computng and writing out LTA {}'.format(ltaFileName))


def save_coregistered_template(imageToImageTransformMatrix, imageToWorldTransformMatrix, savePath, template,
                               templateFileNameBase, templateFileNameExtension):
    # For historical reasons, we applied the estimated transformation to the template; let's do that now
    desiredTemplateImageToWorldTransformMatrix = np.asfortranarray(
        imageToWorldTransformMatrix @ imageToImageTransformMatrix)
    transformedTemplateFileName = os.path.join(savePath,
                                               templateFileNameBase + '_coregistered' + templateFileNameExtension)
    template.write(transformedTemplateFileName, GEMS2Python.KvlTransform(desiredTemplateImageToWorldTransformMatrix))
    return transformedTemplateFileName


if __name__ == '__main__':
    import os

    affineRegistrationMeshCollectionFileName = '/Users/ys/Downloads/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/atlasForAffineRegistration.txt.gz'
    templateFileName = '/Users/ys/Downloads/innolitics_testing/atlas/20Subjects_smoothing2_down2_smoothingForAffine2/template.nii'
    rootPath = '/Users/ys/Downloads/innolitics_testing/buckner40/'
    for f in os.listdir(rootPath):
        fullpath = os.path.join(rootPath, f)
        if os.path.isdir(fullpath):
            imageFileName = os.path.join(fullpath, 'orig.mgz')
            savePath = fullpath
            samseg_registerAtlas(imageFileName=imageFileName,
                                 meshCollectionFileName=affineRegistrationMeshCollectionFileName,
                                 templateFileName=templateFileName,
                                 savePath=savePath,
                                 visualizer=None,
                                 worldToWorldTransformMatrix=None,
                                 InitLTAFile=None,
                                 checkpoint_manager=None)
