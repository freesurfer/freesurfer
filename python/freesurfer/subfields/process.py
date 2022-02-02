import os
import numpy as np
import freesurfer as fs

from freesurfer import samseg
from freesurfer.subfields import utils
from freesurfer.subfields.thalamus import ThalamicNuclei
from freesurfer.subfields.brainstem import BrainstemSubstructures
from freesurfer.subfields.hippocampus import HippoAmygdalaSubfields


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
        fs.fatal(f'Unknown structure type `{structure}`. Available options are: {options}.')
    return model_class


def run_cross_sectional(structure, parameters):
    """
    Run full cross-sectional processing
    """
    total_timer = fs.utils.Timer()

    # Construct model
    model = get_model_class(structure)(**parameters)

    # Preprocess
    timer = fs.utils.Timer()
    print(f'Step 1: preprocessing inputs for {structure} segmentation')
    model.initialize()
    print(f'Preprocessing took {timer.elapsed.seconds} seconds\n')

    # Atlas alignment
    timer = fs.utils.Timer()
    print('Step 2: aligning atlas to reference segmentation')
    model.align_atlas_to_seg()
    print(f'Initial atlas alignment took {timer.elapsed.seconds} seconds\n')

    # Fit to seg
    timer = fs.utils.Timer()
    print('Step 3: fitting mesh to reference segmentation')
    model.prepare_for_seg_fitting()
    model.fit_mesh_to_seg()
    print(f'Initial mesh fitting took {timer.elapsed.seconds} seconds\n')

    # Fit to image
    timer = fs.utils.Timer()
    print('Step 4: fitting mesh to image')
    model.prepare_for_image_fitting()
    model.fit_mesh_to_image()
    print(f'Mesh fitting took {timer.elapsed.seconds} seconds\n')

    # Finalize
    model.extract_segmentation()
    model.postprocess_segmentation()
    model.cleanup()
    print(f'\nSegmentation complete! Process took {total_timer.elapsed.seconds / 60:.2f} minutes')


def run_longitudinal(structure, baseParameters, tpParameters):
    """
    Run full longitudinal processing
    """
    total_timer = fs.utils.Timer()

    # Construct base and timepoint models
    modelClass = get_model_class(structure)
    baseModel = modelClass(**baseParameters)
    tpModels = [modelClass(**params) for params in tpParameters]

    # Preprocess inputs
    timer = fs.utils.Timer()
    print(f'Step 1: preprocessing all inputs for {structure} segmentation')
    baseModel.initialize()
    for tpModel in tpModels:
        tpModel.initialize()
    print(f'Preprocessing took {timer.elapsed.seconds} seconds\n')

    # Atlas alignment
    timer = fs.utils.Timer()
    print('Step 2: aligning atlas to base reference segmentation')
    baseModel.align_atlas_to_seg()
    print(f'Initial base atlas alignment took {timer.elapsed.seconds} seconds\n')

    # Here we compute an affine aligment of the time points to the base based
    # solely on the segmentation. We need to remember the determinant of the
    # transform - we'll need to divide the final volume estimates by it.
    timer = fs.utils.Timer()
    print('Step 3: aligning timepoint segmentations to base')

    # First save the cropped base masks
    mask = baseModel.atlasAlignmentTarget.crop_to_bbox(margin=6)
    mask.data = mask.data.astype('float32') * 255
    baseMaskFile = os.path.join(baseModel.tempDir, 'binaryMaskCropped.mgz')
    mask.write(baseMaskFile)

    # TODOC
    baseProcessedFile = os.path.join(baseModel.tempDir, 'processedImage.mgz')
    baseModel.processedImage.write(baseProcessedFile)

    # Now align each TP to the base
    baseTransforms = []
    for tpModel in tpModels:
        # Save the cropped timepoint masks
        mask = tpModel.atlasAlignmentTarget.crop_to_bbox(margin=6)
        mask.data = mask.data.astype('float32') * 255
        maskFile = os.path.join(tpModel.tempDir, 'binaryMaskCropped.mgz')
        mask.write(maskFile)
        
        # Run the actual registration and load the transform
        transformFile = os.path.join(tpModel.tempDir, 'toBase.lta')
        movedFile = os.path.join(tpModel.tempDir, 'alignedToBase.mgz')
        utils.run(f'mri_robust_register --mov {maskFile} --dst {baseMaskFile} --lta {transformFile} --mapmovhdr {movedFile} --affine --sat 50 -verbose 0')
        baseTransforms.append(fs.LinearTransform.read(transformFile))

        # Resample the inputs in processed base space
        # Again, since we don't have cubic interpolation yet in the python utils, let's just use mri_convert
        correctedImageFile = os.path.join(tpModel.tempDir, 'correctedImage.mgz')
        resampledFile = os.path.join(tpModel.tempDir, 'resampledImage.mgz')
        tpModel.correctedImages[0].write(correctedImageFile)  # TODOC about multi-files
        utils.run(f'mri_convert {correctedImageFile} {resampledFile} -odt float -rl {baseProcessedFile} -rt cubic -at {transformFile}')
        tpModel.processedImage = fs.Volume.read(resampledFile)

        # TODOC
        tpModel.alignedAtlas = baseModel.alignedAtlas

    print(f'Timepoint alignment took {timer.elapsed.seconds} seconds\n')

    # Okay, now we fit the mesh to the base segmentation
    timer = fs.utils.Timer()
    print('Step 4: fitting mesh to base segmentation')
    baseModel.prepare_for_seg_fitting()
    baseModel.fit_mesh_to_seg()
    print(f'Initial mesh fitting took {timer.elapsed.seconds} seconds\n')

    # Global loop: atlas estimation
    timer = fs.utils.Timer()
    print('Step 5: global mesh fitting')

    # TODOC
    baseModel.prepare_for_image_fitting()
    atlasPositions = baseModel.meshCollection.get_mesh(-1).points
    subjectAtlasPositions = baseModel.mesh.points

    # TODOC
    for tpModel in tpModels:
        tpModel.warpedMeshFileName = baseModel.warpedMeshFileName
        tpModel.originalAlphas = baseModel.originalAlphas
        tpModel.prepare_for_image_fitting()
        # TODOC
        tpModel.workingMask = baseModel.workingMask
        tpModel.maskIndices = baseModel.maskIndices
        tpModel.workingImage.data[tpModel.workingMask.data == 0] = 0

    # TODOC
    subjectTPpositions = [tpModel.mesh.points for tpModel in tpModels]

    # This will be the mesh collection we'll use for temporary data
    meshCollection = samseg.gems.KvlMeshCollection()
    meshCollection.read(baseModel.warpedMeshFileName)
    meshCollection.transform(baseModel.transform)
    # meshCollection.K = Katl TODO
    mesh = meshCollection.get_mesh(0)
    mesh.alphas = baseModel.reducedAlphas

    # TODO
    Katl = baseModel.meshStiffness
    Ktp = baseModel.meshStiffness

    # TODOC
    maxGlobalLongIterations = 5  # TODO should be 5!!
    for globalIteration in range(maxGlobalLongIterations):

        print(f'\nGlobal iteration {globalIteration + 1}: estimating subject-specific atlas\n')

        # TODOC
        meshCollection.set_positions(atlasPositions, [subjectAtlasPositions])
        meshSA = meshCollection.get_mesh(0)
        meshCollection.set_positions(atlasPositions, subjectTPpositions)

        # Get optimizer and plug calculator into it
        calculator = samseg.gems.KvlCostAndGradientCalculator(meshCollection, Katl, Ktp, baseModel.transform)
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
                fs.error('minLogLikelihoodTimesPrior is NaN')

            # Check if we need to stop
            if maximalDeformation <= maximalDeformationStopCriterion:
                print('maximalDeformation is too small - stopping')
                break

        # TODO
        # ATH why are we setting K again?
        subjectAtlasPositions = meshSA.points
        baseModel.meshCollection.set_positions(atlasPositions, [subjectAtlasPositions])
        # baseModel.meshCollection.K = K?? TODO
        for t, tpModel in enumerate(tpModels):
            tpModel.meshCollection.set_positions(subjectAtlasPositions, [subjectTPpositions[t]])
            # tpModel.meshCollection.K = K?? TODO
    
        # TODOC
        for t, tpModel in enumerate(tpModels):
            print(f'\nGlobal iteration {globalIteration + 1}: deforming time point {t + 1}\n')

            # TODOC
            idx = 0 if globalIteration < 2 else 1
            tpModel.meshSmoothingSigmas = tpModel.longMeshSmoothingSigmas[idx]
            tpModel.imageSmoothingSigmas = tpModel.longImageSmoothingSigmas[idx]
            tpModel.maxIterations = tpModel.longMaxIterations[idx]

            # TODOC
            tpModel.fit_mesh_to_image()
            # also, check_global_IT iteration init is different when ==1 !

        # Get updated positions
        subjectTPpositions = [tpModel.mesh.points for tpModel in tpModels]

    print(f'Global mesh fitting took {timer.elapsed.seconds} seconds\n')

    # Finalize results
    for tpModel in tpModels:
        tpModel.extract_segmentation()
        
        # TODO here we want to set the header back to where we got it from...

        # Do the subclass-defined postprocessing and cleanup
        tpModel.postprocess_segmentation()
        tpModel.cleanup()

    # Cleanup the base as well
    baseModel.cleanup()
    print(f'\nSegmentation complete! Process took {total_timer.elapsed.seconds / 60:.2f} minutes')
