import freesurfer.gems as gems
import numpy as np
import pytest
import scipy.io

MATLAB_FIXTURE_PATH = 'Testing/matlab_data/integration.mat'
MESHCOLLECTION_PATH = 'Testing/test.txt.gz'
TEST_IMAGE_PATH = 'Testing/test.nii'


@pytest.fixture(scope='session')
def matlab_fixture():
    return scipy.io.loadmat(MATLAB_FIXTURE_PATH)


@pytest.fixture()
def image_fixture():
    return gems.KvlImage(TEST_IMAGE_PATH)


@pytest.fixture()
def cropped_image_fixture():
    # Two different images would be a better test
    return gems.KvlImage(TEST_IMAGE_PATH, TEST_IMAGE_PATH)


@pytest.fixture()
def mesh_collection_fixture():
    collection = gems.KvlMeshCollection()
    collection.read(MESHCOLLECTION_PATH)
    return collection


@pytest.fixture()
def mesh_fixture(mesh_collection_fixture):
    return mesh_collection_fixture.get_mesh(-1)


@pytest.mark.images
@pytest.mark.slowtest
def test_cropped_image_image_bindings(matlab_fixture, cropped_image_fixture):
    np.testing.assert_allclose(matlab_fixture['imageBuffer'], cropped_image_fixture.getImageBuffer())
    expected_transform = matlab_fixture['imageToWorldTransformMatrix']
    np.testing.assert_allclose(expected_transform, cropped_image_fixture.transform_matrix.as_numpy_array)
    assert [189.0, 200.0, 283.0] == cropped_image_fixture.non_cropped_image_size
    assert [0.0, 0.0, 0.0] == cropped_image_fixture.cropping_offset


@pytest.mark.images
@pytest.mark.slowtest
def test_image_bindings(matlab_fixture, image_fixture):
    np.testing.assert_allclose(matlab_fixture['imageBuffer'], image_fixture.getImageBuffer())
    np.testing.assert_allclose(matlab_fixture['imageToWorldTransformMatrix'],
                               image_fixture.transform_matrix.as_numpy_array)


@pytest.mark.images
@pytest.mark.slowtest
def test_smooth_image_buffer(matlab_fixture, image_fixture):
    SIGMA = 1.0
    sigmas = [2 * SIGMA, 3 * SIGMA, 5 * SIGMA]
    input_buffer = image_fixture.getImageBuffer()
    smoothed_buffer = gems.KvlImage.smooth_image_buffer(input_buffer, sigmas)
    assert input_buffer.shape == smoothed_buffer.shape


@pytest.mark.slowtest
def test_mesh_binding(matlab_fixture, image_fixture, mesh_fixture):
    buffer = image_fixture.getImageBuffer()
    np.testing.assert_allclose(matlab_fixture['nodePositions'], mesh_fixture.points)
    np.testing.assert_allclose(matlab_fixture['alphas'], mesh_fixture.alphas)
    rasterized_priors = mesh_fixture.rasterize(buffer.shape)
    assert (189.0, 200.0, 283.0, 17) == rasterized_priors.shape
    np.testing.assert_allclose(matlab_fixture['priors'], rasterized_priors)
    np.testing.assert_allclose(matlab_fixture['nodePositions_after_modification'], mesh_fixture.points + 1)
    mesh_fixture.scale([1 / 10])
    np.testing.assert_allclose(matlab_fixture['nodePositions_after_scaling'], mesh_fixture.points)

@pytest.mark.slowtest
def test_mesh_binding_with_class_name(matlab_fixture, image_fixture, mesh_fixture):
    buffer = image_fixture.getImageBuffer()
    rasterized_priors = mesh_fixture.rasterize(buffer.shape, 0)
    assert (189.0, 200.0, 283.0) == rasterized_priors.shape


@pytest.mark.slowtest
def test_calculator_and_optimzer_bindings(matlab_fixture, image_fixture, mesh_fixture):
    calculator = gems.KvlCostAndGradientCalculator('MutualInformation', [image_fixture], 'Affine')
    mutualInformation_cost, mutualInformation_gradient = calculator.evaluate_mesh_position(mesh_fixture)
    np.testing.assert_almost_equal(matlab_fixture['mutualInformation_cost'][0][0], mutualInformation_cost)
    np.testing.assert_allclose(matlab_fixture['mutualInformation_gradient'], mutualInformation_gradient)

    optimizerType = 'L-BFGS'
    optimizer = gems.KvlOptimizer(optimizerType, mesh_fixture, calculator, {
        'Verbose': 1,
        'MaximalDeformationStopCriterion': 0.1,
        'LineSearchMaximalDeformationIntervalStopCriterion': 0.1,
        'BFGS-MaximumMemoryLength': 12
    })
    minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer()
    np.testing.assert_almost_equal(matlab_fixture['minLogLikelihoodTimesPrior'][0][0], minLogLikelihoodTimesPrior)
    np.testing.assert_almost_equal(matlab_fixture['maximalDeformation'][0][0], maximalDeformation)
