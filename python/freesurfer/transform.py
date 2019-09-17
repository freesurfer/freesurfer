import numpy as np


class LinearTransform:

    def __init__(self, matrix=np.eye(4)):
        self.matrix = matrix

    def transform(self, points):
        pts = np.array(points, copy=False)
        if pts.ndim == 1:
            pts = pts[np.newaxis]
        pts = np.c_[pts, np.ones(pts.shape[0])].T
        return np.dot(self.matrix, pts).T.squeeze()[..., :-1]

    def inverse(self):
        return LinearTransform(np.linalg.inv(self.matrix))
