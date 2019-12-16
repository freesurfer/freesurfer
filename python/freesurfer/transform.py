import os
import numpy as np
from . import bindings


class LinearTransform:
    '''Linear transformer that wraps an affine matrix.'''

    class Type:
        vox = 0
        ras = 1
        physvox = 2

    def __init__(self, matrix=None, source=None, target=None, type=None):
        '''
        Constructs a linear transform from a 4x4 affine matrix. If
        no matrix is provided, the identity is used.
        '''
        self.matrix = matrix if matrix is not None else np.eye(4)
        self.source = source
        self.target = target
        self.type = type

    @staticmethod
    def ensure(x):
        '''TODOC'''
        if isinstance(x, LinearTransform):
            return x
        if isinstance(x, np.ndarray):
            return LinearTransform(x)
        raise ValueError('Cannot convert object of type %s to LinearTransform' % type(unknown).__name__)

    @staticmethod
    def matmul(a, b):
        '''TODOC'''
        a = LinearTransform.ensure(a).matrix
        b = LinearTransform.ensure(b).matrix
        return LinearTransform(np.matmul(a, b))

    def transform(self, points):
        '''Applies the transform matrix to a point or array of points.'''
        pts = np.array(points, copy=False)
        if pts.ndim == 1:
            pts = pts[np.newaxis]
        pts = np.c_[pts, np.ones(pts.shape[0])].T
        return np.dot(self.matrix, pts).T.squeeze()[..., :-1]

    def inverse(self):
        '''Computes the inverse linear transform.'''
        return LinearTransform(np.linalg.inv(self.matrix))

    def write(self, filename):
        bindings.transform.write_lta(self, filename)

    @classmethod
    def read(cls, filename):
        if not os.path.isfile(filename):
            raise ValueError('transform file %s does not exist' % filename)
        return bindings.transform.read_lta(filename)


def LIA(shape=None, voxsize=(1, 1, 1)):
    matrix = np.array([[-1,  0,  0,  0],
                       [ 0,  0,  1,  0],
                       [ 0, -1,  0,  0],
                       [ 0,  0,  0,  1]], dtype=float)
    matrix[:3, :3] *= voxsize
    if shape is not None:
        pcrs = np.append(np.array(shape) / 2, 1)
        matrix[:3, 3] = -np.matmul(matrix, pcrs)[:3]
    return matrix


class Geometry:
    '''TODOC'''

    def __init__(self, shape, voxsize=(1, 1, 1), affine=None):
        self.shape = shape
        self.voxsize = voxsize
        self.affine = affine

    def __eq__(self, other):
        equal = self.shape == other.shape and \
                self.voxsize == other.voxsize and \
                np.array_equal(self.affine, other.affine)
        return equal

    def vox2ras(self):
        '''LinearTransform that maps voxel crs coordinates to ras xyz coordinates.'''
        return LinearTransform(self.affine)

    def ras2vox(self):
        '''LinearTransform that maps ras xyz coordinates to voxel crs coordinates.'''
        return self.vox2ras().inverse()

    def vox2surf(self):
        '''LinearTransform that maps voxel crs coordinates to surface coordinates.'''
        return LinearTransform(LIA(self.shape, self.voxsize))

    def surf2vox(self):
        '''LinearTransform that maps surface coordinates to crs coordinates.'''
        return self.vox2surf().inverse()

    def surf2ras(self):
        '''LinearTransform that maps surface coordinates to ras xyz coordinates.'''
        return LinearTransform.matmul(self.affine, self.surf2vox())

    def ras2surf(self):
        '''LinearTransform that maps ras xyz coordinates to surface coordinates.'''
        return LinearTransform.matmul(self.vox2surf(), self.ras2vox())


class Transformable:
    '''TODOC'''

    def geometry(self):
        raise NotImplementedError('TODOC')

    def vox2ras(self):
        return self.geometry().vox2ras()

    def ras2vox(self):
        return self.geometry().ras2vox()

    def vox2surf(self):
        return self.geometry().vox2surf()

    def surf2vox(self):
        return self.geometry().surf2vox()

    def surf2ras(self):
        return self.geometry().surf2ras()

    def ras2surf(self):
        return self.geometry().ras2surf()

