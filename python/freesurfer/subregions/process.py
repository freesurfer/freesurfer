import os
import numpy as np
import datetime as dt
import surfa as sf

from freesurfer import samseg
from freesurfer.subregions import utils
from freesurfer.subregions.thalamus import ThalamicNuclei
from freesurfer.subregions.brainstem import BrainstemSubstructures
from freesurfer.subregions.hippocampus import HippoAmygdalaSubfields


model_lookup = {
    'thalamus' : ThalamicNuclei,
    'brainstem' : BrainstemSubstructures,
    'hippo-amygdala' : HippoAmygdalaSubfields,
}

structure_names = list(model_lookup.keys())


def get_model_class(structure):
    """
    Get model class from structure name
    """
    model_class = model_lookup.get(structure)
    if model_class is None:
        options = ', '.join(model_lookup.keys())
        sf.system.fatal(f'Unknown structure type `{structure}`. Available options are: {options}.')
    return model_class


def run_cross_sectional(structure, parameters):
    """
    Run full cross-sectional processing
    """
    start_time = dt.datetime.now()

    # Construct model
    model = get_model_class(structure)(**parameters)

    # Preprocess
    last_time = dt.datetime.now()
    print(f'Step 1: preprocessing inputs for {structure} segmentation')
    model.initialize()
    elapsed = dt.datetime.now() - last_time
    print(f'Preprocessing took {elapsed.seconds} seconds\n')

    # Atlas alignment
    last_time = dt.datetime.now()
    print('Step 2: aligning atlas to reference segmentation')
    model.align_atlas_to_seg()
    elapsed = dt.datetime.now() - last_time
    print(f'Initial atlas alignment took {elapsed.seconds} seconds\n')

    # Fit to seg
    last_time = dt.datetime.now()
    print('Step 3: fitting mesh to reference segmentation')
    model.prepare_for_seg_fitting()
    model.fit_mesh_to_seg()
    elapsed = dt.datetime.now() - last_time
    print(f'Initial mesh fitting took {elapsed.seconds} seconds\n')

    # Fit to image
    last_time = dt.datetime.now()
    print('Step 4: fitting mesh to image')
    model.prepare_for_image_fitting()
    model.fit_mesh_to_image()
    elapsed = dt.datetime.now() - last_time
    print(f'Mesh fitting took {elapsed.seconds} seconds\n')

    # Finalize
    model.extract_segmentation()
    model.postprocess_segmentation()
    model.cleanup()
    elapsed = dt.datetime.now() - start_time
    print(f'\nSegmentation complete! Process took {elapsed.seconds / 60:.2f} minutes')


def run_longitudinal(structure, baseParameters, tpParameters):
    """
    Run full longitudinal processing
    """
    start_time = dt.datetime.now()

    # Construct base and timepoint models
    modelClass = get_model_class(structure)
    baseModel = modelClass(**baseParameters)
    tpModels = [modelClass(**params) for params in tpParameters]

    # Preprocess inputs
    last_time = dt.datetime.now()
    print(f'Step 1: preprocessing all inputs for {structure} segmentation')
    baseModel.initialize()
    for tpModel in tpModels:
        tpModel.initialize()
    elapsed = dt.datetime.now() - last_time
    print(f'Preprocessing took {elapsed.seconds} seconds\n')

    # Atlas alignment
    last_time = dt.datetime.now()
    print('Step 2: aligning atlas to base reference segmentation')
    baseModel.align_atlas_to_seg()
    elapsed = dt.datetime.now() - last_time
    print(f'Initial base atlas alignment took {elapsed.seconds} seconds\n')

    # Here we compute an affine aligment of the time points to the base based
    # solely on the segmentation. We need to remember the determinant of the
    # transform - we'll need to divide the final volume estimates by it.
    last_time = dt.datetime.now()
    print('Step 3: aligning timepoint segmentations to base')

    # First save the cropped base masks
    mask = baseModel.atlasAlignmentTarget.crop_to_bbox(margin=6)
    mask.data = mask.data.astype('float32') * 255
    baseMaskFile = os.path.join(baseModel.tempDir, 'binaryMaskCropped.mgz')
    mask.save(baseMaskFile)

    # This is our target resampled image that all timepoints should be resampled to
    baseProcessedFile = os.path.join(baseModel.tempDir, 'processedImage.mgz')
    baseModel.processedImage.save(baseProcessedFile)

    # Now align each TP to the base
    baseTransforms = []
    for tpModel in tpModels:
        # Save the cropped timepoint masks
        mask = tpModel.atlasAlignmentTarget.crop_to_bbox(margin=6)
        mask.data = mask.data.astype('float32') * 255
        maskFile = os.path.join(tpModel.tempDir, 'binaryMaskCropped.mgz')
        mask.save(maskFile)
        
        # Run the actual registration and load the transform
        transformFile = os.path.join(tpModel.tempDir, 'toBase.lta')
        movedFile = os.path.join(tpModel.tempDir, 'alignedToBase.mgz')
        utils.run(f'mri_robust_register --mov {maskFile} --dst {baseMaskFile} --lta {transformFile} --mapmovhdr {movedFile} --affine --sat 50 -verbose 0')
        baseTransforms.append(sf.load_affine(transformFile))

        # Resample the inputs in processed base space
        # Again, since we don't have cubic interpolation yet in the python utils, let's just use mri_convert
        correctedImageFile = os.path.join(tpModel.tempDir, 'correctedImage.mgz')
        resampledFile = os.path.join(tpModel.tempDir, 'resampledImage.mgz')
        tpModel.correctedImages[0].save(correctedImageFile)  # ATH this will need to be adapted for multi-image inputs
        utils.run(f'mri_convert {correctedImageFile} {resampledFile} -odt float -rl {baseProcessedFile} -rt cubic -at {transformFile}')
        tpModel.processedImage = sf.load_volume(resampledFile)

        # Since we're now working in base-space, we can reuse the base-aligned atlas for every timepoint
        tpModel.alignedAtlas = baseModel.alignedAtlas

    elapsed = dt.datetime.now() - last_time
    print(f'Timepoint alignment took {elapsed.seconds} seconds\n')

    # Okay, now we fit the mesh to the base segmentation
    last_time = dt.datetime.now()
    print('Step 4: fitting mesh to base segmentation')
    baseModel.prepare_for_seg_fitting()
    baseModel.fit_mesh_to_seg()
    elapsed = dt.datetime.now() - last_time
    print(f'Initial mesh fitting took {elapsed.seconds} seconds\n')

    # Global loop: atlas estimation
    last_time = dt.datetime.now()
    print('Step 5: global mesh fitting')

    # Prepare the base for image fitting so that we can extract some mesh information
    baseModel.prepare_for_image_fitting()
    atlasPositions = baseModel.meshCollection.get_mesh(-1).points
    subjectAtlasPositions = baseModel.mesh.points

    # Now that we've the mesh to the base segmentation mask, we
    # should use this mesh collection in the timepoint models
    for tpModel in tpModels:
        tpModel.warpedMeshFileName = baseModel.warpedMeshFileName
        tpModel.originalAlphas = baseModel.originalAlphas
        tpModel.prepare_for_image_fitting()

        # We should keep the masking consistent across timepoints
        tpModel.workingMask = baseModel.workingMask
        tpModel.maskIndices = baseModel.maskIndices
        tpModel.workingImage.data[tpModel.workingMask.data == 0] = 0

    # Gather initial subject timepoint mesh positions
    subjectTPpositions = [tpModel.mesh.points for tpModel in tpModels]

    # This will be the mesh collection we'll use for temporary data
    meshCollection = samseg.gems.KvlMeshCollection()
    meshCollection.read(baseModel.warpedMeshFileName)
    meshCollection.transform(baseModel.transform)
    meshCollection.k = baseModel.meshStiffness
    mesh = meshCollection.get_mesh(0)
    mesh.alphas = baseModel.reducedAlphas

    # Start the global iterations
    maxGlobalLongIterations = 5
    for globalIteration in range(maxGlobalLongIterations):

        print(f'\nGlobal iteration {globalIteration + 1}: estimating subject-specific atlas\n')

        # Set temporary data
        meshCollection.set_positions(atlasPositions, [subjectAtlasPositions])
        meshSA = meshCollection.get_mesh(0)
        meshCollection.set_positions(atlasPositions, subjectTPpositions)

        # Get optimizer and plug calculator into it
        calculator = samseg.gems.KvlCostAndGradientCalculator(meshCollection, baseModel.meshStiffness, baseModel.meshStiffness, baseModel.transform)
        maximalDeformationStopCriterion = 1e-10
        optimizationParameters = {
            'Verbose': 0,
            'MaximalDeformationStopCriterion': maximalDeformationStopCriterion,
            'LineSearchMaximalDeformationIntervalStopCriterion': 1e-10,
            'MaximumNumberOfIterations': 1000,
            'BFGS-MaximumMemoryLength': 12
        }
        optimizer = samseg.gems.KvlOptimizer(baseModel.optimizerType, meshSA, calculator, optimizationParameters)

        for positionUpdatingIterationNumber in range(400):

            # Calculate a good step
            minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_samseg()

            # Log optimization information
            iteration_info = [
                f'GlobalIter: {globalIteration + 1}',
                f'Iter: {positionUpdatingIterationNumber + 1:03d}',
                f'MaxDef: {maximalDeformation:.4f}',
                f'MinLLxP: {minLogLikelihoodTimesPrior:.4f}',
            ]
            print('  '.join(iteration_info))

            if np.isnan(minLogLikelihoodTimesPrior):
                print('error: minLogLikelihoodTimesPrior is NaN')

            # Check if we need to stop
            if maximalDeformation <= maximalDeformationStopCriterion:
                print('maximalDeformation is too small - stopping')
                break

        # Update positions
        subjectAtlasPositions = meshSA.points
        baseModel.meshCollection.set_positions(atlasPositions, [subjectAtlasPositions])
        for t, tpModel in enumerate(tpModels):
            tpModel.meshCollection.set_positions(subjectAtlasPositions, [subjectTPpositions[t]])
    
        # Now let's fit each timepoint
        for t, tpModel in enumerate(tpModels):
            print(f'\nGlobal iteration {globalIteration + 1}: deforming time point {t + 1}\n')

            # Update the multi-resolution settings
            idx = 0 if globalIteration < 2 else 1
            tpModel.meshSmoothingSigmas = tpModel.longMeshSmoothingSigmas[idx]
            tpModel.imageSmoothingSigmas = tpModel.longImageSmoothingSigmas[idx]
            tpModel.maxIterations = tpModel.longMaxIterations[idx]

            # Do the fitting stage
            tpModel.fit_mesh_to_image()

        # Get updated positions
        subjectTPpositions = [tpModel.mesh.points for tpModel in tpModels]

    elapsed = dt.datetime.now() - last_time
    print(f'Global mesh fitting took {elapsed.seconds} seconds\n')

    # Finalize results
    for t, tpModel in enumerate(tpModels):
        tpModel.extract_segmentation()
        # Let's transform (just the header) the output segmentations back to original timepoint space
        trf = baseTransforms[t]
        tpModel.discreteLabels.geom.vox2world = trf.inv() @ tpModel.discreteLabels.geom.vox2world
        # Also, scale the volumes by the determinant of the transform
        det = np.linalg.det(trf.matrix[:3, :3])
        tpModel.volumes = {key: vol / det for key, vol in tpModel.volumes.items()}
        # Do the subclass-defined postprocessing and cleanup
        tpModel.postprocess_segmentation()
        tpModel.cleanup()

    # Cleanup the base as well
    baseModel.cleanup()
    elapsed = dt.datetime.now() - start_time
    print(f'\nSegmentation complete! Process took {elapsed.seconds / 60:.2f} minutes')
