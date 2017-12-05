import logging
import os

from easydict import EasyDict

logger = logging.getLogger(__name__)


def update_recipe_with_calculated_paths(recipe, avg_data_dir=None):
    if avg_data_dir is None:
        avg_data_dir = os.environ.get('SAMSEG_DATA_DIR')
    recipe.save_path = find_or_create_save_path(recipe)
    recipe.mesh_collection_file_name = determine_mesh_collection_file_name(avg_data_dir)
    recipe.compression_lookup_table_file_name = determine_compression_lookup_table_file_name(avg_data_dir)
    recipe.show_segmentation_figures = recipe.exvivo
    recipe.show_registration_figures = False
    recipe.template_file_name = determine_template_file_name(avg_data_dir)


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
