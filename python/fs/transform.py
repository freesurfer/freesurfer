import numpy as np


class LinearTransform:
    '''Linear transformer that wraps an affine matrix.'''

    def __init__(self, matrix=None):
        '''Constructs a linear transform from a 4x4 affine matrix. If
        no matrix is provided, the identity is used.'''
        self.matrix = matrix if matrix else np.eye(4)

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
