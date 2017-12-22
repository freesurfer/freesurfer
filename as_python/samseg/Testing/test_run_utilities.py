import numpy as np

from as_python.samseg.run_utilities import find_or_create_save_path, determine_compression_lookup_table_file_name, \
    determine_optimization_options, specify_model, Specification

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


class TestFindOrCreateSavePath:
    def setup(self):
        self.made_at = None
        self.okay = None

    def test_find_or_create_save_path(self):
        recipe = Specification({'output': 'write_here'})

        def makedirs(path, exist_ok=False):
            self.made_at = path
            self.okay = exist_ok
            return path

        save_path = find_or_create_save_path(recipe, makedirs=makedirs)
        assert 'write_here' == save_path
        assert 'write_here' == self.made_at
        assert self.okay


def test_determine_optimization_options():
    for verbose in [True, False]:
        actual = determine_optimization_options(verbose, avg_data_dir='yogi')
        actual_multi_resolution_specification = actual.multiResolutionSpecification
        assert 2 == len(actual_multi_resolution_specification)
        for index, spec in enumerate(actual_multi_resolution_specification):
            assert 'yogi/atlas_level{0}.txt.gz'.format(index + 1) == spec.atlasFileName
            assert 100 == spec.maximumNumberOfIterations
            assert spec.estimateBiasField
            assert 2.0 == spec.targetDownsampledVoxelSpacing or index != 0
            assert 1.0 == spec.targetDownsampledVoxelSpacing or index != 1
        assert 20 == actual.maximumNumberOfDeformationIterations
        assert 1e-4 == actual.absoluteCostPerVoxelDecreaseStopCriterion
        assert verbose == actual.verbose
        assert 0.001 == actual.maximalDeformationStopCriterion
        assert actual.maximalDeformationStopCriterion == \
               actual.lineSearchMaximalDeformationIntervalStopCriterion
        assert 0.0 == actual.maximalDeformationAppliedStopCriterion
        assert 12 == actual.BFGSMaximumMemoryLength


def test_specify_model():
    shared_gmm_parameters = 'hello'
    FreeSurferLabels = range(5)
    colors = [1, 2, 3]
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
