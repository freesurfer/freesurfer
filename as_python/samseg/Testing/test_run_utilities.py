import pytest
from easydict import EasyDict

from as_python.samseg.run_utilities import find_or_create_save_path, determine_mesh_collection_file_name, \
    determine_compression_lookup_table_file_name, determine_transformed_template_filename, determine_template_file_name, \
    exvivo_shared_gmm_parameters, standard_shared_gmm_parameters

EXPECTED_NAMES = {name for name in [
    'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus',
    'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm',
    'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities', 'Putamen',
    'Thalamus','Pallidum',
    ]}

class TestRunUtilities:
    def setup(self):
        self.made_at = None
        self.okay = None

    def test_determine_compression_lookup_table_file_name(self):
        actual = determine_compression_lookup_table_file_name('panda')
        expected = 'panda/namedCompressionLookupTable.txt'
        assert expected == actual

    def test_determine_mesh_collection_file_name(self):
        actual = determine_mesh_collection_file_name('grizzly')
        expected = 'grizzly/CurrentMeshCollection30New.txt.gz'
        assert expected == actual

    def test_determine_mesh_collection_file_name(self):
        actual = determine_mesh_collection_file_name('black')
        expected = 'black/CurrentMeshCollection30New.txt.gz'
        assert expected == actual

    def test_determine_template_file_name(self):
        actual = determine_template_file_name('teddy')
        expected = 'teddy/mni305_masked_autoCropped.mgz'
        assert expected == actual

    def test_determine_transformed_template_filename(self):
        actual = determine_transformed_template_filename('polar')
        expected = 'polar/mni305_masked_autoCropped_coregistered.mgz'
        assert expected == actual

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

    def test_exvivo_shared_gmm_parameters(self):
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

    def test_standard_shared_gmm_parameters(self):
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

def name_set(gmm):
    return {name
            for item in gmm
            for name in item.search_strings}