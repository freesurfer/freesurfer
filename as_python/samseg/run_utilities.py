import logging
import os

import numpy as np
from easydict import EasyDict

logger = logging.getLogger(__name__)

def find_avg_data_dir(avg_data_dir):
    if avg_data_dir is None:
        return os.environ.get('SAMSEG_DATA_DIR')
    return avg_data_dir

def update_recipe_with_calculated_paths(recipe, avg_data_dir=None):
    avg_data_dir = find_avg_data_dir(avg_data_dir)
    recipe = EasyDict(recipe)
    recipe.save_path = find_or_create_save_path(recipe)
    recipe.mesh_collection_file_name = determine_mesh_collection_file_name(avg_data_dir)
    recipe.compression_lookup_table_file_name = determine_compression_lookup_table_file_name(avg_data_dir)
    recipe.show_segmentation_figures = recipe.exvivo
    recipe.show_registration_figures = False
    recipe.template_file_name = determine_template_file_name(avg_data_dir)
    recipe.avg_data_dir = avg_data_dir
    return recipe


def determine_compression_lookup_table_file_name(avg_data_dir):
    return '{0}/compressionLookupTable.txt'.format(avg_data_dir)


def determine_mesh_collection_file_name(avg_data_dir):
    return '{0}/CurrentMeshCollection30New.txt.gz'.format(avg_data_dir)


def determine_template_file_name(avg_data_dir):
    return '{0}/mni305_masked_autoCropped.mgz'.format(avg_data_dir)


def determine_transformed_template_filename(save_path):
    return '{0}/mni305_masked_autoCropped_coregistered.mgz'.format(save_path)


def find_or_create_save_path(recipe, makedirs=os.makedirs):
    save_path = recipe.output
    makedirs(save_path, exist_ok=True)
    return save_path


def determine_shared_gmm_parameters(exvivo):
    if exvivo:
        return exvivo_shared_gmm_parameters()
    else:
        return standard_shared_gmm_parameters()


def exvivo_shared_gmm_parameters():
    return [
        #   % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
        #   sharedGMMParameters = struct;
        #   sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background and what is normally CSF
        #   sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
        #   sharedGMMParameters( 1 ).numberOfComponents = 1;
        EasyDict({
            'mergedName': 'Unknown',
            'searchStrings': ['Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus'],
            'numberOfComponents': 1,
        }),
        #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
        #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
        #   sharedGMMParameters( 2 ).numberOfComponents = 1;
        EasyDict({
            'mergedName': 'Global WM',
            'searchStrings': ['White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm'],
            'numberOfComponents': 1,
        }),
        #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
        #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen' };
        #   sharedGMMParameters( 3 ).numberOfComponents = 1;
        EasyDict({
            'mergedName': 'Global GM',
            'searchStrings': ['Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities',
                               'Putamen'],
            'numberOfComponents': 1,
        }),
        #   sharedGMMParameters( 4 ).mergedName = 'Thalamus'; % Thalamus
        #   sharedGMMParameters( 4 ).searchStrings = { 'Thalamus' };
        #   sharedGMMParameters( 4 ).numberOfComponents = 1;
        EasyDict({
            'mergedName': 'Thalamus',
            'searchStrings': ['Thalamus'],
            'numberOfComponents': 1,
        }),
        #   sharedGMMParameters( 5 ).mergedName = 'Pallidum'; % Pallidum
        #   sharedGMMParameters( 5 ).searchStrings = { 'Pallidum' };
        #   sharedGMMParameters( 5 ).numberOfComponents = 1;
        EasyDict({
            'mergedName': 'Pallidum',
            'searchStrings': ['Pallidum'],
            'numberOfComponents': 1,
        }),
    ]


def standard_shared_gmm_parameters():
    return [
        #   % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
        #   sharedGMMParameters = struct;
        #   sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background
        #   sharedGMMParameters( 1 ).searchStrings = { 'Unknown'};
        #   sharedGMMParameters( 1 ).numberOfComponents = 3;
        EasyDict({
            'mergedName': 'Unknown',
            'searchStrings': ['Unknown'],
            'numberOfComponents': 2,
        }),
        #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
        #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
        #   sharedGMMParameters( 2 ).numberOfComponents = 2;
        EasyDict({
            'mergedName': 'Global WM',
            'searchStrings': ['White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm'],
            'numberOfComponents': 2,
        }),
        #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
        #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities' };
        #   sharedGMMParameters( 3 ).numberOfComponents = 3;
        EasyDict({
            'mergedName': 'Global GM',
            'searchStrings': ['Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities'],
            'numberOfComponents': 2,
        }),
        #   sharedGMMParameters( 4 ).mergedName = 'Global CSF'; % CSF
        #   sharedGMMParameters( 4 ).searchStrings = { 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
        #   sharedGMMParameters( 4 ).numberOfComponents = 3;
        EasyDict({
            'mergedName': 'Global CSF',
            'searchStrings': ['Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus'],
            'numberOfComponents': 3,
        }),
        #   sharedGMMParameters( 5 ).mergedName = 'Thalamus'; % Thalamus
        #   sharedGMMParameters( 5 ).searchStrings = { 'Thalamus' };
        #   sharedGMMParameters( 5 ).numberOfComponents = 2;
        EasyDict({
            'mergedName': 'Thalamus',
            'searchStrings': ['Thalamus'],
            'numberOfComponents': 2,
        }),
        #   sharedGMMParameters( 6 ).mergedName = 'Pallidum'; % Pallidum
        #   sharedGMMParameters( 6 ).searchStrings = { 'Pallidum' };
        #   sharedGMMParameters( 6 ).numberOfComponents = 2;
        EasyDict({
            'mergedName': 'Pallidum',
            'searchStrings': ['Pallidum'],
            'numberOfComponents': 2,
        }),
        #   sharedGMMParameters( 7 ).mergedName = 'Putamen'; % Putamen
        #   sharedGMMParameters( 7 ).searchStrings = { 'Putamen' };
        #   sharedGMMParameters( 7 ).numberOfComponents = 2;
        EasyDict({
            'mergedName': 'Putamen',
            'searchStrings': ['Putamen'],
            'numberOfComponents': 2,
        }),
    ]


def determine_optimization_options(verbose=False, avg_data_dir=None):
    avg_data_dir = find_avg_data_dir(avg_data_dir)
    return EasyDict({
        'multiResolutionSpecification': [
            # % Set various optimization options
            # optimizationOptions = struct;
            # optimizationOptions.multiResolutionSpecification = struct;
            # optimizationOptions.multiResolutionSpecification(1).atlasFileName = fullfile(samsegDataDir, 'atlas_level1.txt.gz');
            # optimizationOptions.multiResolutionSpecification( 1 ).meshSmoothingSigma = 2.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 1 ).estimateBiasField = true;
            {
                'atlasFileName': os.path.join(avg_data_dir, 'atlas_level1.txt.gz'),
                'meshSmoothingSigma': 2.0,
                'targetDownsampledVoxelSpacing': 2.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True,
            },
            # optimizationOptions.multiResolutionSpecification( 2 ).atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
            # optimizationOptions.multiResolutionSpecification( 2 ).meshSmoothingSigma = 0.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 2 ).estimateBiasField = true; % Switching this off will use the bias field estimated
            #                                                                                 % at lower resolution(s)
            {
                'atlasFileName': os.path.join(avg_data_dir, 'atlas_level2.txt.gz'),
                'meshSmoothingSigma': 0.0,
                'targetDownsampledVoxelSpacing': 1.0,
                'maximumNumberOfIterations': 100,
                'estimateBiasField': True,
            },
        ],
        # optimizationOptions.maximumNumberOfDeformationIterations = 20;
        'maximumNumberOfDeformationIterations': 20,
        # optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
        'absoluteCostPerVoxelDecreaseStopCriterion': 1e-4,
        # optimizationOptions.verbose = 0;
        'verbose': verbose,
        # optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
        'maximalDeformationStopCriterion': 0.001,
        # optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
        'lineSearchMaximalDeformationIntervalStopCriterion': 0.001,
        # % optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
        'relativeCostDecreaseStopCriterion': 1e-6,
        # optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
        'maximalDeformationAppliedStopCriterion': 0.0,
        # optimizationOptions.BFGSMaximumMemoryLength = 12;
        'BFGSMaximumMemoryLength': 12,
    })


def specify_model(FreeSurferLabels, noBrainMasking, useDiagonalCovarianceMatrices, shared_gmm_parameters, names, colors, avg_data_dir=None):
    avg_data_dir = find_avg_data_dir(avg_data_dir)
    return EasyDict({
        'FreeSurferLabels': FreeSurferLabels,
        # modelSpecifications.atlasFileName = fullfile( samsegDataDir, 'atlas_level2.txt.gz' );
        'atlasFileName': os.path.join(avg_data_dir, 'atlas_level2.txt.gz'),
        # modelSpecifications.FreeSurferLabels = FreeSurferLabels;
        # modelSpecifications.names = names;
       'names': names,
        # modelSpecifications.colors = colors;
        'colors': colors,
        # modelSpecifications.sharedGMMParameters = sharedGMMParameters;
        'sharedGMMParameters': shared_gmm_parameters,
        # modelSpecifications.useDiagonalCovarianceMatrices = useDiagonalCovarianceMatrices;
        'useDiagonalCovarianceMatrices': useDiagonalCovarianceMatrices,
        # modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel
        'brainMaskingSmoothingSigma': 3.0,
        # if noBrainMasking
        #   modelSpecifications.brainMaskingThreshold = -Inf;
        # else
        #   modelSpecifications.brainMaskingThreshold = 0.01;
        'brainMaskingThreshold': -np.inf if noBrainMasking else 0.01,
        # end
        # modelSpecifications.K = 0.1; % Stiffness of the mesh
        'K': 0.1,
        # modelSpecifications.biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing
        'biasFieldSmoothingKernelSize': 50,
    })


def use_standard_affine_registration_atlas(avg_data_dir):
    #   affineRegistrationMeshCollectionFileName = sprintf( '%s/SPM12_6classes_30x30x30_meshCollection.txt.gz', AvgDataDir );
    #   affineRegistrationTemplateFileName = sprintf( '%s/SPM12_6classes_30x30x30_template.nii', AvgDataDir );
    return EasyDict({
        'mesh_collection_file_name':
            '{0}/atlasForAffineRegistration.txt.gz'.format(avg_data_dir),
        'template_file_name':
            '{0}/template.nii'.format(avg_data_dir),
    })
