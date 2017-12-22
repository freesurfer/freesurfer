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
            for name in item.searchStrings}

def test_determine_compression_lookup_table_file_name():
    actual = determine_compression_lookup_table_file_name('panda')
    expected = 'panda/compressionLookupTable.txt'
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
        assert 1 == item.numberOfComponents
        assert len(item.searchStrings) > 0
    assert 5 == len(actual)
    assert 'Unknown' == actual[0].mergedName
    assert 'Global WM' == actual[1].mergedName
    assert 'Global GM' == actual[2].mergedName
    assert 'Thalamus' == actual[3].mergedName
    assert 'Pallidum' == actual[4].mergedName
    assert ['Cortex',
            'Caudate',
            'Hippocampus',
            'Amygdala',
            'Accumbens',
            'hypointensities',
            'Putamen'] == actual[2].searchStrings
    assert EXPECTED_NAMES == name_set(actual)

def test_standard_shared_gmm_parameters():
    actual = standard_shared_gmm_parameters()
    assert 7 == len(actual)
    assert 'Unknown' == actual[0].mergedName
    assert 'Global WM' == actual[1].mergedName
    assert 'Global GM' == actual[2].mergedName
    assert 'Global CSF' == actual[3].mergedName
    assert 'Thalamus' == actual[4].mergedName
    assert 'Pallidum' == actual[5].mergedName
    assert 'Putamen' == actual[6].mergedName
    assert ['Cortex',
            'Caudate',
            'Hippocampus',
            'Amygdala',
            'Accumbens',
            'hypointensities',
            ] == actual[2].searchStrings
    for index, item in enumerate(actual):
        actual_number_of_components = item.numberOfComponents
        assert actual_number_of_components == 3 or index != 3
        assert actual_number_of_components == 2 or index == 3
        assert len(item.searchStrings) > 0
    assert EXPECTED_NAMES == name_set(actual)

def test_determine_shared_gmm_parameters():
    assert 7 == len(determine_shared_gmm_parameters(False))
    assert 5 == len(determine_shared_gmm_parameters(True))

def test_determine_optimization_options():
    for verbose in [True, False]:
        actual = determine_optimization_options(verbose, avg_data_dir='yogi')
        actual_multi_resolution_specification = actual.multiResolutionSpecification
        assert 2 == len(actual_multi_resolution_specification)
        for index, spec in enumerate(actual_multi_resolution_specification):
            assert 'yogi/atlas_level{0}.txt.gz'.format(index + 1) == spec.atlasFileName
            assert 100 == spec.maximumNumberOfIterations
            assert spec.estimateBiasField
            assert 2.0 == spec.meshSmoothingSigma or index != 0
            assert 0.0 == spec.meshSmoothingSigma or index != 1
            assert 2.0 == spec.targetDownsampledVoxelSpacing or index != 0
            assert 1.0 == spec.targetDownsampledVoxelSpacing or index != 1
        assert 20 == actual.maximumNumberOfDeformationIterations
        assert 1e-4 == actual.absoluteCostPerVoxelDecreaseStopCriterion
        assert verbose == actual.verbose
        assert 0.001 == actual.maximalDeformationStopCriterion
        assert actual.maximalDeformationStopCriterion == \
               actual.lineSearchMaximalDeformationIntervalStopCriterion
        assert 1e-6 == actual.relativeCostDecreaseStopCriterion
        assert 0.0 == actual.maximalDeformationAppliedStopCriterion
        assert 12 == actual.BFGSMaximumMemoryLength

def test_specify_model():
    shared_gmm_parameters = standard_shared_gmm_parameters()
    FreeSurferLabels = range(5)
    colors = [1,2,3]
    names = ['cat', 'dog']
    for noBrainMasking in [True, False]:
        for useDiagonalCovarianceMatrices in [True, False]:
            actual = specify_model(
                FreeSurferLabels,
                noBrainMasking,
                useDiagonalCovarianceMatrices,
                shared_gmm_parameters,
                names,
                colors,
                avg_data_dir='booboo')
            assert 'booboo/atlas_level2.txt.gz' == actual.atlasFileName
            assert useDiagonalCovarianceMatrices == actual.useDiagonalCovarianceMatrices
            assert shared_gmm_parameters == actual.sharedGMMParameters
            assert FreeSurferLabels == actual.FreeSurferLabels
            assert 3.0 == actual.brainMaskingSmoothingSigma
            assert 0.1 == actual.K
            assert 50.0 == actual.biasFieldSmoothingKernelSize
            assert 0.01 == actual.brainMaskingThreshold or noBrainMasking
            assert -np.inf == actual.brainMaskingThreshold or not noBrainMasking
            assert colors == actual.colors
            assert names == actual.names
