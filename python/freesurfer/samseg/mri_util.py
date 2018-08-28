import gzip

import nibabel as nib
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
    image = nib.load(filename)
    header = image.header
    volsz = header['dims']
    if header['goodRASFlag']:
        M = header.get_affine()
    else:
        M = None
    return [None, M, header['mrparms'], np.array(volsz)]


def construct_affine(mat, offset):
    M = np.zeros([4, 4])
    for row in range(3):
        for column in range(3):
            M[row, column] = mat[row, column]
        M[row, 3] = offset[row]
    M[3, 3] = 1.0
    return M
