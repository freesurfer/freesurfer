import GEMS2Python
import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal

from as_python.samseg.kvl import transform_product, voxel_spacing_of_transform, calculate_down_sampling_factors


def test_transform_product():
    matrix_a = np.array([
        [10.0, 2, 3, 1],
        [5, 20, 7, 2],
        [11, 13, 30, 3],
        [0, 0, 0, 1],
    ], dtype=np.double, order='F')
    matrix_b = np.array([
        [100.0, 17, 19, 10],
        [23, 200, 29, 20],
        [31, 37, 300, 30],
        [0, 0, 0, 1],
    ], dtype=np.double, order='F')

    transform_a = GEMS2Python.KvlTransform(matrix_a)
    transform_b = GEMS2Python.KvlTransform(matrix_b)

    matrix_ab = np.dot(matrix_a, matrix_b)
    transform_ab = transform_product(transform_a, transform_b)

    actual_ab = transform_ab.as_numpy_array
    assert_array_equal(matrix_ab, actual_ab)


def test_voxel_spacing_of_transform():
    test_matrix = GEMS2Python.KvlTransform(np.array([
        [200, 300, 500, 999],
        [20, 30, 50, 999],
        [2, 3, 5, 999],
        [0, 0, 0, 1],
    ], dtype=np.double, order='F'))
    actual_spacing = voxel_spacing_of_transform(test_matrix)
    assert 3 == len(actual_spacing)
    actual_spacing_squared = [v*v for v in actual_spacing]
    expected_spacing_squared = [(1*1 + 10*10 + 100*100)*v for v in [2*2, 3*3, 5*5]]
    assert_array_almost_equal(expected_spacing_squared, actual_spacing_squared)


def test_calculate_down_sampling_factors():
    test_matrix = GEMS2Python.KvlTransform(np.array([
        [10, 0, 0, 0],
        [0, 20, 0, 0],
        [0, 0, 30, 0],
        [0, 0, 0, 1],
    ], dtype=np.double, order='F'))
    actual_small = calculate_down_sampling_factors(test_matrix, 10)
    assert_array_equal([1, 1, 1], actual_small)
    actual_big = calculate_down_sampling_factors(test_matrix, 135)
    assert_array_equal([14, 7, 4], actual_big)
