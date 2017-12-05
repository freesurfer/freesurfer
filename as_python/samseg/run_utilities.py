import logging
import os
import numpy as np

from easydict import EasyDict

logger = logging.getLogger(__name__)


def update_recipe_with_calculated_paths(recipe, avg_data_dir=None):
    recipe = EasyDict(recipe)
    if avg_data_dir is None:
        avg_data_dir = os.environ.get('SAMSEG_DATA_DIR')
    recipe.save_path = find_or_create_save_path(recipe)
    recipe.mesh_collection_file_name = determine_mesh_collection_file_name(avg_data_dir)
    recipe.compression_lookup_table_file_name = determine_compression_lookup_table_file_name(avg_data_dir)
    recipe.show_segmentation_figures = recipe.exvivo
    recipe.show_registration_figures = False
    recipe.template_file_name = determine_template_file_name(avg_data_dir)
    recipe.avg_data_dir = avg_data_dir
    return recipe


def determine_compression_lookup_table_file_name(avg_data_dir):
    return '{0}/namedCompressionLookupTable.txt'.format(avg_data_dir)


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
            'merged_name': 'Unknown',
            'search_strings': ['Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus'],
            'number_of_components': 1,
        }),
        #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
        #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
        #   sharedGMMParameters( 2 ).numberOfComponents = 1;
        EasyDict({
            'merged_name': 'Global WM',
            'search_strings': ['White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm'],
            'number_of_components': 1,
        }),
        #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
        #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen' };
        #   sharedGMMParameters( 3 ).numberOfComponents = 1;
        EasyDict({
            'merged_name': 'Global GM',
            'search_strings': ['Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities',
                               'Putamen'],
            'number_of_components': 1,
        }),
        #   sharedGMMParameters( 4 ).mergedName = 'Thalamus'; % Thalamus
        #   sharedGMMParameters( 4 ).searchStrings = { 'Thalamus' };
        #   sharedGMMParameters( 4 ).numberOfComponents = 1;
        EasyDict({
            'merged_name': 'Thalamus',
            'search_strings': ['Thalamus'],
            'number_of_components': 1,
        }),
        #   sharedGMMParameters( 5 ).mergedName = 'Pallidum'; % Pallidum
        #   sharedGMMParameters( 5 ).searchStrings = { 'Pallidum' };
        #   sharedGMMParameters( 5 ).numberOfComponents = 1;
        EasyDict({
            'merged_name': 'Pallidum',
            'search_strings': ['Pallidum'],
            'number_of_components': 1,
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
            'merged_name': 'Unknown',
            'search_strings': ['Unknown'],
            'number_of_components': 2,
        }),
        #   sharedGMMParameters( 2 ).mergedName = 'Global WM'; % WM
        #   sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
        #   sharedGMMParameters( 2 ).numberOfComponents = 2;
        EasyDict({
            'merged_name': 'Global WM',
            'search_strings': ['White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm'],
            'number_of_components': 2,
        }),
        #   sharedGMMParameters( 3 ).mergedName = 'Global GM'; % GM
        #   sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities' };
        #   sharedGMMParameters( 3 ).numberOfComponents = 3;
        EasyDict({
            'merged_name': 'Global GM',
            'search_strings': ['Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities'],
            'number_of_components': 2,
        }),
        #   sharedGMMParameters( 4 ).mergedName = 'Global CSF'; % CSF
        #   sharedGMMParameters( 4 ).searchStrings = { 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
        #   sharedGMMParameters( 4 ).numberOfComponents = 3;
        EasyDict({
            'merged_name': 'Global CSF',
            'search_strings': ['Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus'],
            'number_of_components': 3,
        }),
        #   sharedGMMParameters( 5 ).mergedName = 'Thalamus'; % Thalamus
        #   sharedGMMParameters( 5 ).searchStrings = { 'Thalamus' };
        #   sharedGMMParameters( 5 ).numberOfComponents = 2;
        EasyDict({
            'merged_name': 'Thalamus',
            'search_strings': ['Thalamus'],
            'number_of_components': 2,
        }),
        #   sharedGMMParameters( 6 ).mergedName = 'Pallidum'; % Pallidum
        #   sharedGMMParameters( 6 ).searchStrings = { 'Pallidum' };
        #   sharedGMMParameters( 6 ).numberOfComponents = 2;
        EasyDict({
            'merged_name': 'Pallidum',
            'search_strings': ['Pallidum'],
            'number_of_components': 2,
        }),
        #   sharedGMMParameters( 7 ).mergedName = 'Putamen'; % Putamen
        #   sharedGMMParameters( 7 ).searchStrings = { 'Putamen' };
        #   sharedGMMParameters( 7 ).numberOfComponents = 2;
        EasyDict({
            'merged_name': 'Putamen',
            'search_strings': ['Putamen'],
            'number_of_components': 2,
        }),
    ]


def determine_optimization_options(verbose=False):
    return EasyDict({
        'multi_resolution_specification': [
            # % Set various optimization options
            # optimizationOptions = struct;
            # optimizationOptions.multiResolutionSpecification = struct;
            # optimizationOptions.multiResolutionSpecification( 1 ).meshSmoothingSigma = 2.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 1 ).targetDownsampledVoxelSpacing = 2.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 1 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 1 ).estimateBiasField = true;
            {
                'mesh_smoothing_sigma': 2.0,
                'target_downsampling_voxel_spacing': 2.0,
                'maximum_number_of_iterations': 100,
                'estimate_bias_field': True,
            },
            # optimizationOptions.multiResolutionSpecification( 2 ).meshSmoothingSigma = 0.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 2 ).targetDownsampledVoxelSpacing = 1.0; % In mm
            # optimizationOptions.multiResolutionSpecification( 2 ).maximumNumberOfIterations = 100;
            # optimizationOptions.multiResolutionSpecification( 2 ).estimateBiasField = true; % Switching this off will use the bias field estimated
            #                                                                                 % at lower resolution(s)
            {
                'mesh_smoothing_sigma': 0.0,
                'target_downsampling_voxel_spacing': 1.0,
                'maximum_number_of_iterations': 100,
                'estimate_bias_field': True,
            },
        ],
        # optimizationOptions.maximumNumberOfDeformationIterations = 20;
        'maximum_number_of_deformation_iterations': 20,
        # optimizationOptions.absoluteCostPerVoxelDecreaseStopCriterion = 1e-4;
        'absolute_cost_per_voxel_decreases_stop_criterion': 1e-4,
        # optimizationOptions.verbose = 0;
        'verbose': verbose,
        # optimizationOptions.maximalDeformationStopCriterion = 0.001; % Measured in pixels
        'maximal_deformation_stop_criterion': 0.001,
        # optimizationOptions.lineSearchMaximalDeformationIntervalStopCriterion = optimizationOptions.maximalDeformationStopCriterion; % Idem
        'line_search_maximal_deformation_interval_stop_criterion': 0.001,
        # % optimizationOptions.relativeCostDecreaseStopCriterion = 1e-6;
        'relative_cost_decrease_stop_criterion': 1e-6,
        # optimizationOptions.maximalDeformationAppliedStopCriterion = 0.0;
        'maximal_deformation_applied_stop_criterion': 0.0,
        # optimizationOptions.BFGSMaximumMemoryLength = 12;
        'bfgs_maximum_memory_length': 12,
    })


def specify_model(exvivo, missing_structures, shared_gmm_parameters):
    return EasyDict({
        # exvivo = cmdargs.exvivo;
        # % Set various model specifications
        # modelSpecifications = struct;
        # modelSpecifications.missingStructureSearchStrings = missingStructureSearchStrings;
        'missing_structures': missing_structures,
        # modelSpecifications.sharedGMMParameters = sharedGMMParameters;
        'shared_gmm_parameters': shared_gmm_parameters,
        # modelSpecifications.useDiagonalCovarianceMatrices = false;
        'use_diagonal_covariance_matrices': exvivo,
        # modelSpecifications.brainMaskingSmoothingSigma = 3; % sqrt of the variance of a Gaussian blurring kernel
        'brain_masking_smoothing_sigma': 3.0,
        # modelSpecifications.brainMaskingThreshold = 0.01;
        'brain_masking_threshold': -np.inf if exvivo else 0.01,
        # modelSpecifications.K = 0.1; % Stiffness of the mesh
        'k': 0.1,
        # modelSpecifications.biasFieldSmoothingKernelSize = 50.0;  % Distance in mm of sinc function center to first zero crossing
        'bias_field_smoothing_kernel_size': 50,
        # if exvivo
        #   modelSpecifications.brainMaskingThreshold = -Inf; % Disable brain masking
        #   modelSpecifications.useDiagonalCovarianceMatrices = true;
        # end
    })


def use_standard_affine_registration_atlas(avg_data_dir):
    #   affineRegistrationMeshCollectionFileName = sprintf( '%s/SPM12_6classes_30x30x30_meshCollection.txt.gz', AvgDataDir );
    #   affineRegistrationTemplateFileName = sprintf( '%s/SPM12_6classes_30x30x30_template.nii', AvgDataDir );
    return EasyDict({
        'mesh_collection_file_name':
            '{0}/SPM12_6classes_30x30x30_meshCollection.txt.gz'.format(avg_data_dir),
        'template_file_name':
            '{0}/SPM12_6classes_30x30x30_template.nii'.format(avg_data_dir),
    })
