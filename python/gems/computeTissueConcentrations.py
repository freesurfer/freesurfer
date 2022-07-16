#
# Compute tissue concentrations maps by altering the posterior maps obtained with Samseg
# with the Jacobian determinant of the deformation.
#
# Command example:
# --subjects-dir /path/to/my/samseg/runs
# --output-dir /path/to/output
# --template /usr/local/freesurfer/7.2.0/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/template.nii
# --mesh-collection-file /usr/local/freesurfer/7.2.0/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/atlas_level2.txt.gz
# --structures "Left-Cerebral-Cortex" "Right-Cerebral-Cortex"
# --downsampling-factors 1 2 3
# --save-figs
# --use-block-smoothing
#
# Add the flag --longitudinal for runs obtained with SamsegLongitudinal
#
# Note that all runs must have been executed with the --history flag
# and with either the command run_samseg or run_samseg_long (not samseg or samseg_long!)
#

import os
import numpy as np
import argparse
import surfa as sf
from scipy.ndimage import map_coordinates
from freesurfer import samseg

parser = argparse.ArgumentParser()
parser.add_argument('--subjects-dir', help='Directory with saved SAMSEG runs with --history flag.', required=True)
parser.add_argument('--output-dir', help='Output directory.', required=True)
parser.add_argument('--template', help='Template file name.', required=True)
parser.add_argument('--mesh-collection-file', help='Mesh collection file name.', required=True)
parser.add_argument('--structures', nargs='+', help='Structure(s) of interest.', required=True)
parser.add_argument('--downsampling-factors', nargs='+', type=int, help='Downsampling factor(s).', required=True)
parser.add_argument('--use-block-smoothing', action='store_true', default=False, help='Use block smoothing instead of Gaussian smoothing.')
parser.add_argument('--longitudinal', action='store_true', default=False, help='Longitudinal SAMSEG runs.')
parser.add_argument('--show-figs', action='store_true', default=False, help='Show figures during run.')
parser.add_argument('--save-figs', action='store_true', default=False, help='Save averages.')

args = parser.parse_args()

# Make output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

if args.show_figs:
    visualizer = samseg.initVisualizer(True, True)
else:
    visualizer = samseg.initVisualizer(False, False)

if args.save_figs:
    import nibabel as nib

# We need an init of the probabilistic segmentation class
# to call instance methods
atlas = samseg.ProbabilisticAtlas()

subjectList = [pathname for pathname in os.listdir(args.subjects_dir) if os.path.isdir(os.path.join(args.subjects_dir, pathname))]
subjectList.sort()
numberOfSubjects = len(subjectList)

print('Number of subjects: ' + str(len(subjectList)))

# Load template beforehand as this is common for all the subjects
template = nib.load(args.template)
templateBuffer = template.get_fdata()
templateTransformMatrix = template.affine

# Same for mesh collection with reference mesh
mesh_collection_file = args.mesh_collection_file
meshCollection = samseg.gems.KvlMeshCollection()
meshCollection.read(mesh_collection_file)
referenceMesh = meshCollection.reference_mesh
referencePositions = referenceMesh.points.copy()
numberOfNodes = referenceMesh.point_count

# Avoid boundary issues with mesh rasterizor for tetrahedra that have faces that are truly perfectly
# aligned with the Z-plane
referencePositions = referencePositions + 0.01 * np.random.random([numberOfNodes, 3])
referenceMesh.points = referencePositions

# Concentration maps are what we are interest in
# I'm using lists here to avoid having to have many if else when dealing with longitudinal runs
concentrationMaps = []

for structure in args.structures:

    print("Working on structure: " + str(structure))

    for subjectNumber, subjectDir in enumerate(subjectList):

        print("Working on subject: " + str(subjectNumber + 1))

        concentrationMapsSubject = []

        # Read the structure posterior
        if not args.longitudinal:
            posteriorPaths = [os.path.join(args.subjects_dir, subjectDir, "posteriors", str(structure) + ".mgz")]
            historyPaths = [os.path.join(args.subjects_dir, subjectDir, "history.p")]
        else:
            timePoints = os.listdir(os.path.join(args.subjects_dir, subjectDir))
            timePoints.sort()
            timePoints.remove('base')
            timePoints.remove('latentAtlases')
            timePoints.remove('history.p')
            history_latent = np.load(os.path.join(args.subjects_dir, subjectDir, "history.p"), allow_pickle=True)
            latent_deformation = history_latent['latentDeformationEvolution'][-1]
            posteriorPaths = [os.path.join(args.subjects_dir, subjectDir, tpDir, "posteriors", str(structure) + ".mgz") for tpDir in timePoints]
            historyPaths = [os.path.join(args.subjects_dir, subjectDir, tpDir, "history.p") for tpDir in timePoints]

        for posteriorPath, historyPath in zip(posteriorPaths, historyPaths):
            segmentationMap = nib.load(posteriorPath).get_fdata()

            # Load history file with all the information
            history = np.load(historyPath, allow_pickle=True)
            model_specifications = history['input']['modelSpecifications']
            transform_matrix = history['transform']
            transform = samseg.gems.KvlTransform(samseg.requireNumpyArray(transform_matrix))
            if not args.longitudinal:
                deformation = history['historyWithinEachMultiResolutionLevel'][1]['deformation']
            else:
                deformation = history['historyWithinEachMultiResolutionLevel'][0]['deformation'] + latent_deformation
            mesh = atlas.getMesh(mesh_collection_file,
                                 transform,
                                 K=model_specifications.K,
                                 initialDeformation=deformation,
                                 initialDeformationMeshCollectionFileName=mesh_collection_file)
            positions = mesh.points.copy()
            # The image is cropped as well so the voxel coordinates
            # do not exactly match with the original image,
            # i.e., there's a shift. Let's undo that.
            cropping = history['cropping']
            positions += [slc.start for slc in cropping]

            # Alter the tissue probability map with the Jacobian determinant of the deformation --
            # resulting in the so - called "tissue concentration" maps
            mesh.points = positions.copy()
            expansionMap = mesh.draw_jacobian_determinant(segmentationMap.shape)

            #
            segmentationMap = segmentationMap * expansionMap


            if args.show_figs:
                visualizer.show(images=segmentationMap)

            #
            minPositions = np.min(positions, axis=0)
            maxPositions = np.max(positions, axis=0)
            rescaledPositions = (positions - minPositions) / (maxPositions - minPositions)
            referenceMesh.alphas = rescaledPositions
            tmp = referenceMesh.rasterize(templateBuffer.shape, -1)

            # Unvisited voxels are marked by zeroes in all three coordinates. Make sure those only
            # occur around the image boundaries (and not in the middle, which would/could destroy
            # any voxel-based models we might want to use
            validMask = np.sum(tmp != 0, 3) != 0
            if args.show_figs:
                visualizer.show(images=validMask)
            validMaskInside = validMask[1:-1, 1:-1, 1:-1]
            assert np.alltrue(validMaskInside == 1)

            # Now resample - we're (ab)using the fact here that unvisited voxels are mapped to coordinate
            # (1,1,1) in template space, which is a corner where values should be automatically zero in
            # subject space (corner of affinely aligned template box)
            densePositionsX = (tmp[:, :, :, 0] / 65535) * (maxPositions[0] - minPositions[0]) + minPositions[0]
            densePositionsY = (tmp[:, :, :, 1] / 65535) * (maxPositions[1] - minPositions[1]) + minPositions[1]
            densePositionsZ = (tmp[:, :, :, 2] / 65535) * (maxPositions[2] - minPositions[2]) + minPositions[2]

            resampledImageBuffer = map_coordinates(segmentationMap, (densePositionsX, densePositionsY, densePositionsZ),
                                                   order=1, cval=0)

            if args.show_figs:
                visualizer.show(images=resampledImageBuffer)

            # Store results
            concentrationMapsSubject.append(resampledImageBuffer)

        concentrationMaps.append(concentrationMapsSubject)

    if args.show_figs:
        visualizer.show(images=templateBuffer)
        averageConcentrationMap = np.mean(np.vstack(concentrationMaps), axis=0)
        visualizer.show(images=averageConcentrationMap)
    if args.save_figs:
        averageConcentrationMap = np.mean(np.vstack(concentrationMaps), axis=0)
        img = nib.Nifti1Image(averageConcentrationMap, np.eye(4))
        nib.save(img, os.path.join(args.output_dir, str(structure) + "_averageConcentrationMap"))
        img = nib.Nifti1Image(templateBuffer, np.eye(4))
        nib.save(img, os.path.join(args.output_dir, "template"))

    # Also try (smoothed) downsampling
    for downSamplingFactor in args.downsampling_factors:
        print("Computing downsampled maps with downSamplingFactor: " + str(downSamplingFactor))
        size = templateBuffer.shape
        downsampledConcentrationMaps = []

        if not args.use_block_smoothing:
            # Use Gaussian smoothing
            if downSamplingFactor == 1:
                smoothingSigmas = [0.0, 0.0, 0.0]
            else:
                # Variance chosen to approximately match normalized binomial filter
                # (1/4, 1/2, 1/4) for downsampling factor of 2
                smoothingSigmas = (downSamplingFactor / 2 / np.sqrt(2 * np.log(2)) * np.ones([1, 3])).tolist()[0]


            for subjectNumber, subjectDir in enumerate(subjectList):

                downsampledConcentrationMapsSubject = []

                tmps = [np.array(concentrationMaps[subjectNumber][tp], dtype=np.float32, order='F')
                        for tp in range(len(concentrationMaps[subjectNumber]))]

                for tmp in tmps:
                    smoothedConcentrationMap = samseg.gems.KvlImage.smooth_image_buffer(tmp, smoothingSigmas)
                    downsampledConcentrationMap = smoothedConcentrationMap[::downSamplingFactor, ::downSamplingFactor, ::downSamplingFactor]
                    downsampledConcentrationMapsSubject.append(downsampledConcentrationMap)

                downsampledConcentrationMaps.append(downsampledConcentrationMapsSubject)

            smoothedTemplateBuffer = samseg.gems.KvlImage.smooth_image_buffer(templateBuffer, smoothingSigmas)
            downsampledTemplateBuffer = smoothedTemplateBuffer[::downSamplingFactor, ::downSamplingFactor, ::downSamplingFactor]
        else:
            # Use Block smoothing
            for subjectNumber, subjectDir in enumerate(subjectList):

                downsampledConcentrationMapsSubject = []

                tmps = [concentrationMaps[subjectNumber][tp] for tp in range(len(concentrationMaps[subjectNumber]))]
                for tmp in tmps:
                    currentSize = list(tmp.shape)
                    for dimension in range(3):
                        # Reshape into 2-D, do the work in the first dimension, and shape into N-D
                        tmp = np.reshape(tmp, [currentSize[0], -1])
                        downsampledSizeInThisDimension = int(np.ceil(currentSize[0] / downSamplingFactor))
                        # No need to use sparse identity matrix here and scipy.sparse matrices are tricky to use in combination with numpy
                        filterConcentration = np.kron(np.eye(downsampledSizeInThisDimension), np.ones([1, downSamplingFactor]) / downSamplingFactor)
                        filterConcentration = filterConcentration[:, 0:currentSize[0]]
                        tmp = filterConcentration @ tmp
                        currentSize[0] = downsampledSizeInThisDimension
                        tmp = np.reshape(tmp, currentSize)
                        # Shift dimension
                        currentSize = [currentSize[1], currentSize[2], currentSize[0]]
                        tmp = np.moveaxis(tmp, 0, 2)

                    downsampledConcentrationMapsSubject.append(tmp)

                downsampledConcentrationMaps.append(downsampledConcentrationMapsSubject)

            # Same for template
            tmp = templateBuffer
            currentSize = list(tmp.shape)
            for dimension in range(3):
                # Reshape into 2-D, do the work in the first dimension, and shape into N-D
                tmp = np.reshape(tmp, [currentSize[0], -1])
                downsampledSizeInThisDimension = int(np.ceil(currentSize[0] / downSamplingFactor))
                # No need to use sparse identity matrix as sizes are still manageable and scipy.sparse matrices are tricky to use with numpy
                filterConcentration = np.kron(np.eye(downsampledSizeInThisDimension),
                                              np.ones([1, downSamplingFactor]) / downSamplingFactor)
                filterConcentration = filterConcentration[:, 0:currentSize[0]]
                tmp = filterConcentration @ tmp
                currentSize[0] = downsampledSizeInThisDimension
                tmp = np.reshape(tmp, currentSize)
                # Shift dimension
                currentSize = [currentSize[1], currentSize[2], currentSize[0]]
                tmp = np.moveaxis(tmp, 0, 2)

            downsampledTemplateBuffer = tmp

        if args.show_figs:
            averageDownsampledConcentrationMap = np.mean(np.vstack(downsampledConcentrationMaps), axis=0)
            visualizer.show(images=averageDownsampledConcentrationMap)
            visualizer.show(images=downsampledTemplateBuffer)
        if args.save_figs:
            averageDownsampledConcentrationMap = np.mean(np.vstack(downsampledConcentrationMaps), axis=0)

            img = nib.Nifti1Image(averageDownsampledConcentrationMap, np.eye(4))
            nib.save(img, os.path.join(args.output_dir, str(structure) + "_DS_" +
                                       str(downSamplingFactor) + "_averageDownsampledConcentrationMap"))
            img = nib.Nifti1Image(downsampledTemplateBuffer, np.eye(4))
            nib.save(img, os.path.join(args.output_dir, "DS_" +
                                       str(downSamplingFactor) + "_downsampledTemplate"))

        # Now save things
        if not args.use_block_smoothing:
            fileEnding = ".npz"
        else:
            fileEnding = "_blockSmoothed.npz"
        fileName = os.path.join(args.output_dir, str(structure) + "_down_" + str(downSamplingFactor) + fileEnding)

        # If cross-sectional save everything as a numpy array
        if not args.longitudinal:
            downsampledConcentrationMaps = np.vstack(downsampledConcentrationMaps)

        np.savez(fileName, downsampledTemplateBuffer=downsampledTemplateBuffer, downsampledConcentrationMaps=downsampledConcentrationMaps)


# @Stefano:
# Code for using --save-mesh instead of --save-history
# There is something off however, probably due to different saved positions(?) -- is this a bug of --save-mesh?
# Load SAMSEG original mesh
#meshCollection = samseg.gems.KvlMeshCollection()
#meshCollection.read(os.path.join(args.subjects_dir, subjectDir, "mesh.txt.gz"))
#referenceMesh = meshCollection.reference_mesh
#referencePositions = referenceMesh.points.copy()
#numberOfNodes = referenceMesh.point_count
#mesh = meshCollection.get_mesh(0) # Index 0 has the deformed mesh
#positions = mesh.points.copy()

# Avoid boundary issues with mesh rasterizor for tetrahedra that have faces that are truly perfectly
# aligned with the Z-plane
#referencePositions = referencePositions + 0.01 * np.random.random([numberOfNodes, 3])
#referenceMesh.points = referencePositions

# Apply affine transform to the deformed node position
#tmp = np.ones([numberOfNodes, 4])
#tmp[:, 0:3] = positions.copy()
#tmp = (affineTransformMatrix @ tmp.T ).T
#positions = tmp[:, 0:3]

# Alter the tissue probability map with the Jacobian determinant of the deformation --
# resulting in the so - called "tissue concentration" maps
#mesh.points = positions
#expansionMap = mesh.draw_jacobian_determinant(segmentationMap.shape)

# Undo expansion due to global head size
#expansionMap = expansionMap / np.linalg.det(affineTransformMatrix[0:3, 0:3])
#segmentationMap = segmentationMap * expansionMap