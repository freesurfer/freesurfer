import logging
import numpy as np
import gems


logger = logging.getLogger(__name__)


# todo: note: none of this stuff is used

def transform_product(a, b):
    aa = a.as_numpy_array
    bb = b.as_numpy_array
    ab = np.dot(aa, bb)
    ab_f = np.asfortranarray(ab)
    return gems.KvlTransform(ab_f)


def voxel_spacing_of_transform(t):
    m = t.as_numpy_array
    mm = m[0:3, 0:3]
    mm2 = mm * mm
    ss = mm2.sum(axis=0)
    return np.sqrt(ss)


def calculate_down_sampling_factors(transform, target_voxel_spacing):
    # Figure out how much to downsample (depends on voxel size)
    voxel_spacing = voxel_spacing_of_transform(transform)
    tentative_spacing = np.rint(target_voxel_spacing / voxel_spacing)
    limited_spacing = np.maximum([1, 1, 1], tentative_spacing)
    logger.debug('limited_spacing = %s', str(limited_spacing))
    return [int(factor) for factor in limited_spacing]


def create_transform_from_2d_list(list_2d):
    return gems.KvlTransform(np.array(list_2d, dtype=np.double, order='F'))


def create_translation_transform(offsets):
    return create_transform_from_2d_list([
        [1, 0, 0, offsets[0]],
        [0, 1, 0, offsets[1]],
        [0, 0, 1, offsets[2]],
        [0, 0, 0, 1],
    ])
