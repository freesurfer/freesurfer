import numpy as np

from . import bindings


class Geometry:
    def __init__(self, shape, voxsize=(1, 1, 1), affine=None):
        self.shape = shape
        self.voxsize = voxsize
        self.affine = affine


class TransformType:
    vox = 0
    ras = 1
    physvox = 2


class LinearTransform:
    '''Linear transformer that wraps an affine matrix.'''

    def __init__(self, matrix=None, source=None, target=None, type=None):
        '''
        Constructs a linear transform from a 4x4 affine matrix. If
        no matrix is provided, the identity is used.
        '''
        self.matrix = matrix if matrix is not None else np.eye(4)
        self.source = source
        self.target = target
        self.type = type

    def __mul__(self, xform):
        xform = xform.matrix if isinstance(xform, LinearTransform) else xform
        return LinearTransform(np.matmul(self.matrix, xform))  # TODO require

    def __rmul__(self, xform):
        xform = xform.matrix if isinstance(xform, LinearTransform) else xform
        return LinearTransform(np.matmul(xform, self.matrix))

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

    @classmethod
    def read(cls, filename):
        return bindings.read_lta(filename)
