import pytest
import numpy as np
import scipy.io
import GEMS2Python

MATLAB_FIXTURE_PATH = 'Testing/matlab_data/integration.mat'
MESHCOLLECTION_PATH = 'Testing/test.txt.gz'
TEST_IMAGE_PATH = 'Testing/test.nii'

@pytest.fixture(scope='session')
def matlab_fixture():
    return scipy.io.loadmat(MATLAB_FIXTURE_PATH)


@pytest.fixture()
def image_fixture():
    return GEMS2Python.KvlImage(TEST_IMAGE_PATH)


@pytest.fixture()
def mesh_collection_fixture():
    collection = GEMS2Python.KvlMeshCollection()
    collection.read(MESHCOLLECTION_PATH)
    return collection


@pytest.fixture()
def mesh_fixture(mesh_collection_fixture):
    return mesh_collection_fixture.get_mesh(-1)


def test_image_bindings(matlab_fixture, image_fixture):
    np.testing.assert_allclose(matlab_fixture['imageBuffer'], image_fixture.getImageBuffer())
    np.testing.assert_allclose(matlab_fixture['imageToWorldTransformMatrix'], image_fixture.getTransform().as_numpy_array)


def test_mesh_binding(matlab_fixture, image_fixture, mesh_fixture):
    buffer = image_fixture.getImageBuffer()
    np.testing.assert_allclose(matlab_fixture['nodePositions'], mesh_fixture.points)
    np.testing.assert_allclose(matlab_fixture['alphas'], mesh_fixture.alphas)
    rasterized_priors = mesh_fixture.rasterize(buffer.shape)
    np.testing.assert_allclose(matlab_fixture['priors'], rasterized_priors)
    np.testing.assert_allclose(matlab_fixture['nodePositions_after_modification'], mesh_fixture.points + 1)
    mesh_fixture.scale([1/10])
    np.testing.assert_allclose(matlab_fixture['nodePositions_after_scaling'], mesh_fixture.points)


def test_calculator_and_optimzer_bindings(matlab_fixture, image_fixture, mesh_fixture):
    calculator = GEMS2Python.KvlCostAndGradientCalculator( 'MutualInformation', [image_fixture], 'Affine')
    mutualInformation_cost, mutualInformation_gradient = calculator.evaluate_mesh_position(mesh_fixture)
    np.testing.assert_almost_equal(matlab_fixture['mutualInformation_cost'][0][0], mutualInformation_cost)
    np.testing.assert_allclose(matlab_fixture['mutualInformation_gradient'], mutualInformation_gradient)

    optimizerType = 'L-BFGS'
    optimizer = GEMS2Python.KvlOptimizer( optimizerType, mesh_fixture, calculator, {
        'Verbose': 1,
        'MaximalDeformationStopCriterion': 0.1,
        'LineSearchMaximalDeformationIntervalStopCriterion': 0.1,
        'BFGS-MaximumMemoryLength': 12
    })
    minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer()
    np.testing.assert_almost_equal(matlab_fixture['minLogLikelihoodTimesPrior'][0][0], minLogLikelihoodTimesPrior)
    np.testing.assert_almost_equal(matlab_fixture['maximalDeformation'][0][0], maximalDeformation)


