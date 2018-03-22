import gzip
import os
import struct
from functools import reduce
from operator import mul

import numpy as np

FILE_OPENERS = {
    '.mgh': open,
    '.mgz': gzip.open,
}

MRI_UCHAR = 0
MRI_INT = 1
MRI_LONG = 2
MRI_FLOAT = 3
MRI_SHORT = 4
MRI_BITMAP = 5

BYTES_PER_VOXEL = {
    MRI_FLOAT: 4,
    MRI_UCHAR: 1,
    MRI_SHORT: 2,
    MRI_INT: 4,
}


def load_mgh_header(filename):
    with open_mri_file(filename) as file:
        def read_int():
            return struct.unpack('>i', file.read(4))[0]

        def read_int_list(count):
            return [read_int() for _ in range(count)]

        def read_float():
            return struct.unpack('>f', file.read(4))[0]

        def read_float_list(count):
            return [read_float() for _ in range(count)]

        def read_short():
            return struct.unpack('>h', file.read(2))[0]

        vol = None

        v = read_int()
        volsz = read_int_list(4)
        file_type = read_int()
        dof = read_int()

        # read transform, if present
        USED_SPACE_SIZE = (3 * 4 + 4 * 3 * 4)  # spacefor ras transform
        UNUSED_SPACE_SIZE = 256
        unused_space_size = UNUSED_SPACE_SIZE - 2
        ras_good_flag = read_short()
        if ras_good_flag:
            delta = np.array(read_float_list(3))
            Mdc = np.array([read_float_list(3) for _ in range(3)]).T
            Pxyz_c = np.array(read_float_list(3))
            unused_space_size -= USED_SPACE_SIZE
            M = construct_transform(delta, Mdc, Pxyz_c, volsz)
        else:
            M = None

        # skip to start of voxels
        file.seek(unused_space_size, 1)  # from current position

        # skip voxels
        number_of_voxels = reduce(mul, volsz, 1)
        bytes_per_voxel = BYTES_PER_VOXEL[file_type]
        file.seek(number_of_voxels * bytes_per_voxel, 1)  # from current position

        mr_parms = [read_float() for _ in range(4)]
        return [vol, M, mr_parms, volsz]


def construct_transform(delta, Mdc, Pxyz_c, volsz):
    # D = diag(delta);
    D = np.diag(delta)
    # Pcrs_c = [ndim1 / 2 ndim2 / 2 ndim3 / 2]'
    ndims = volsz[:3]
    Pcrs_c = np.array([0.5 * ndim for ndim in ndims]).T  # % Should this be kept?
    # Pxyz_0 = Pxyz_c - Mdc * D * Pcrs_c;
    MdcD = np.matmul(Mdc, D)
    Pxyz_0 = Pxyz_c - np.matmul(MdcD, Pcrs_c)
    # M = [Mdc * D Pxyz_0;
    # ...
    # 0
    # 0
    # 0
    # 1];
    M = np.zeros([4, 4])
    for row in range(3):
        for column in range(3):
            M[row, column] = MdcD[row, column]
        M[row, 3] = Pxyz_0[row]
    M[3, 3] = 1.0
    return M


def open_mri_file(filename):
    _, extension = os.path.splitext(filename)
    opener = FILE_OPENERS.get(extension.lower())
    if opener:
        return opener(filename, 'rb')
    else:
        raise Exception('Unexpected file name extension {}'.format(extension))
