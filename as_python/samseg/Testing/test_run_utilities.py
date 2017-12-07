from easydict import EasyDict
import numpy as np

from as_python.samseg.run_utilities import find_or_create_save_path, determine_mesh_collection_file_name, \
    determine_compression_lookup_table_file_name, determine_transformed_template_filename, determine_template_file_name, \
    exvivo_shared_gmm_parameters, standard_shared_gmm_parameters, determine_optimization_options, specify_model, \
    determine_shared_gmm_parameters, update_recipe_with_calculated_paths, use_standard_affine_registration_atlas

EXPECTED_NAMES = {name for name in [
    'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus',
    'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm',
    'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen',
    'Thalamus', 'Pallidum',
]}

def name_set(gmm):
    return {name
            for item in gmm
            for name in item.search_strings}

def test_determine_compression_lookup_table_file_name():
    actual = determine_compression_lookup_table_file_name('panda')
    expected = 'panda/namedCompressionLookupTable.txt'
    assert expected == actual

def test_determine_mesh_collection_file_name():
    actual = determine_mesh_collection_file_name('grizzly')
    expected = 'grizzly/CurrentMeshCollection30New.txt.gz'
    assert expected == actual

def test_determine_mesh_collection_file_name():
    actual = determine_mesh_collection_file_name('black')
    expected = 'black/CurrentMeshCollection30New.txt.gz'
    assert expected == actual

def test_determine_template_file_name():
    actual = determine_template_file_name('teddy')
    expected = 'teddy/mni305_masked_autoCropped.mgz'
    assert expected == actual

def test_determine_transformed_template_filename():
    actual = determine_transformed_template_filename('polar')
    expected = 'polar/mni305_masked_autoCropped_coregistered.mgz'
    assert expected == actual

def test_update_recipe_with_calculated_paths():
    recipe = EasyDict({'output':'write_here', 'exvivo': True})
    actual = update_recipe_with_calculated_paths(recipe, 'bears')
    assert 'bears' == actual.avg_data_dir

def test_use_standard_affine_registration_atlas():
    actual = use_standard_affine_registration_atlas('affine')
    assert 'affine/atlasForAffineRegistration.txt.gz' == \
           actual.mesh_collection_file_name
    assert 'affine/template.nii' == \
           actual.template_file_name

class TestFindOrCreateSavePath:
    def setup(self):
        self.made_at = None
        self.okay = None


    def test_find_or_create_save_path(self):
        recipe = EasyDict({'output': 'write_here'})

        def makedirs(path, exist_ok=False):
            self.made_at = path
            self.okay = exist_ok
            return path

        save_path = find_or_create_save_path(recipe, makedirs=makedirs)
        assert 'write_here' == save_path
        assert 'write_here' == self.made_at
        assert self.okay

def test_exvivo_shared_gmm_parameters():
    actual = exvivo_shared_gmm_parameters()
    for item in actual:
        assert 1 == item.number_of_components
        assert len(item.search_strings) > 0
    assert 5 == len(actual)
    assert 'Unknown' == actual[0].merged_name
    assert 'Global WM' == actual[1].merged_name
    assert 'Global GM' == actual[2].merged_name
    assert 'Thalamus' == actual[3].merged_name
    assert 'Pallidum' == actual[4].merged_name
    assert ['Cortex',
            'Caudate',
            'Hippocampus',
            'Amygdala',
            'Accumbens',
            'hypointensities',
            'Putamen'] == actual[2].search_strings
    assert EXPECTED_NAMES == name_set(actual)

def test_standard_shared_gmm_parameters():
    actual = standard_shared_gmm_parameters()
    assert 7 == len(actual)
    assert 'Unknown' == actual[0].merged_name
    assert 'Global WM' == actual[1].merged_name
    assert 'Global GM' == actual[2].merged_name
    assert 'Global CSF' == actual[3].merged_name
    assert 'Thalamus' == actual[4].merged_name
    assert 'Pallidum' == actual[5].merged_name
    assert 'Putamen' == actual[6].merged_name
    assert ['Cortex',
            'Caudate',
            'Hippocampus',
            'Amygdala',
            'Accumbens',
            'hypointensities',
            ] == actual[2].search_strings
    for index, item in enumerate(actual):
        actual_number_of_components = item.number_of_components
        assert actual_number_of_components == 3 or index != 3
        assert actual_number_of_components == 2 or index == 3
        assert len(item.search_strings) > 0
    assert EXPECTED_NAMES == name_set(actual)

def test_determine_shared_gmm_parameters():
    assert 7 == len(determine_shared_gmm_parameters(False))
    assert 5 == len(determine_shared_gmm_parameters(True))

def test_determine_optimization_options():
    for verbose in [True, False]:
        actual = determine_optimization_options(verbose)
        actual_multi_resolution_specification = actual.multi_resolution_specification
        assert 2 == len(actual_multi_resolution_specification)
        for index, spec in enumerate(actual_multi_resolution_specification):
            assert 100 == spec.maximum_number_of_iterations
            assert spec.estimate_bias_field
            assert 2.0 == spec.mesh_smoothing_sigma or index != 0
            assert 0.0 == spec.mesh_smoothing_sigma or index != 1
            assert 2.0 == spec.target_downsampling_voxel_spacing or index != 0
            assert 1.0 == spec.target_downsampling_voxel_spacing or index != 1
        assert 20 == actual.maximum_number_of_deformation_iterations
        assert 1e-4 == actual.absolute_cost_per_voxel_decreases_stop_criterion
        assert verbose == actual.verbose
        assert 0.001 == actual.maximal_deformation_stop_criterion
        assert actual.maximal_deformation_stop_criterion == \
               actual.line_search_maximal_deformation_interval_stop_criterion
        assert 1e-6 == actual.relative_cost_decrease_stop_criterion
        assert 0.0 == actual.maximal_deformation_applied_stop_criterion
        assert 12 == actual.bfgs_maximum_memory_length

def test_specify_model():
    shared_gmm_parameters = standard_shared_gmm_parameters()
    missing_structures = ['Cortex', 'Caudate']
    for exvivo in [True, False]:
        actual = specify_model(exvivo, missing_structures, shared_gmm_parameters)
        assert missing_structures == actual.missing_structures
        assert shared_gmm_parameters == actual.shared_gmm_parameters
        assert 3.0 == actual.brain_masking_smoothing_sigma
        assert 0.1 == actual.k
        assert 50.0 == actual.bias_field_smoothing_kernel_size
        assert exvivo == actual.use_diagonal_covariance_matrices
        assert 0.01 == actual.brain_masking_threshold or exvivo
        assert -np.inf == actual.brain_masking_threshold or not exvivo
