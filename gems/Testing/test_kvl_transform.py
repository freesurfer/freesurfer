import freesurfer.gems as gems
import numpy as np

class TestKvlTransform:
    def test_numpy_create(self):
        initial_transform = np.array([
            [1, 2, 3, 4],
            [11, 22, 33, 44],
            [111, 222, 333, 444],
            [0, 0, 0, 1],
        ], dtype=np.double, order='F')
        kvl_transform = gems.KvlTransform(initial_transform)
        returned_transform = kvl_transform.as_numpy_array
        for row in range(4):
            for column in range(4):
                expected = initial_transform[row][column]
                actual = returned_transform[row][column]
                assert expected == actual
