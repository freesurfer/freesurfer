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

        # make sure matrix is 4x4
        if self.matrix.shape == (3, 4):
            tmp = np.eye(4)
            tmp[:3, :] = self.matrix
            self.matrix = tmp

        if self.matrix.shape != (4, 4):
            raise ValueError('linear transform must be a 4x4 matrix')

    @staticmethod
    def ensure(unknown):
        '''Ensures that the unknown input is (or gets converted to) a LinearTransform'''
        if isinstance(unknown, LinearTransform):
            return unknown
        if isinstance(unknown, np.ndarray):
            return LinearTransform(unknown)
        raise ValueError('Cannot convert object of type %s to LinearTransform' % type(unknown).__name__)

    @staticmethod
    def matmul(a, b):
        '''
        Matrix multiplication to join two transforms. Input can be numpy
        arrays or LinearTransforms.
        '''
        a = LinearTransform.ensure(a).matrix
        b = LinearTransform.ensure(b).matrix
        return LinearTransform(np.matmul(a, b))

    def transform(self, points):
        '''Applies the transform matrix to a point or array of points.'''
        pts = np.array(points, copy=False)
        if pts.ndim == 1:
            pts = pts[np.newaxis]
        pts = np.c_[pts, np.ones(pts.shape[0])].T
        return np.ascontiguousarray(np.dot(self.matrix, pts).T.squeeze()[..., :-1])

    def inverse(self):
        '''Computes the inverse linear transform.'''
        return LinearTransform(np.linalg.inv(self.matrix), source=self.target, target=self.source, type=self.type)

    def as_ras(self):
        '''Converts affine matrix to a RAS to RAS transform.'''
        if self.type is None:
            raise ValueError('transform must have type if converting between spaces')
        # if already the desired type, return self
        if self.type == LinearTransform.Type.ras:
            return self
        # include source/target RAS information
        matrix = self.target.affine @ self.matrix @ np.linalg.inv(self.source.affine)
        return LinearTransform(matrix, source=self.source, target=self.target, type=LinearTransform.Type.ras)

    def as_vox(self):
        '''Converts affine matrix to a VOX to VOX transform.'''
        if self.type is None:
            raise ValueError('transform must have type if converting between spaces')
        # if already the desired type, return self
        if self.type == LinearTransform.Type.vox:
            return self
        # exclude source/target RAS information
        matrix = np.linalg.inv(self.target.affine) @ self.matrix @ self.source.affine
        return LinearTransform(matrix, source=self.source, target=self.target, type=LinearTransform.Type.vox)

    def write(self, filename):
        '''Writes the transform to an LTA file.'''
        bindings.transform.write_lta(self, filename)

    @classmethod
    def read(cls, filename):
        '''Reads a transform from an LTA file.'''
        if not os.path.isfile(filename):
            raise ValueError('transform file %s does not exist' % filename)
        return bindings.transform.read_lta(filename)


def LIA(shape=None, voxsize=(1, 1, 1)):
    '''Affine matrix representing LIA voxel coordinate relationship (also known as RSP).'''
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
    '''Representation of volume geometry in RAS space.'''

    def __init__(self, shape, voxsize=(1, 1, 1), affine=None):
        self.shape = shape
        self.voxsize = voxsize
        self.affine = affine

    def __eq__(self, other):
        '''Check for equality.'''
        equal = self.shape == other.shape and \
                self.voxsize == other.voxsize and \
                np.array_equal(self.affine, other.affine)
        return equal

    @staticmethod
    def is_equal(a, b, thresh=1e-3, require_affine=True):
        '''Compare geometries within some threshold.'''

        differ = lambda a, b: not np.allclose(a, b, rtol=0.0, atol=thresh)

        if differ(a.shape, b.shape):
            return False
        if differ(a.voxsize, b.voxsize):
            return False

        if a.affine is not None and b.affine is not None:
            if differ(a.affine, b.affine):
                print(np.abs(np.array(a.affine) - np.array(b.affine)))
                return False
        elif require_affine:
            return False

        return True

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
    '''
    Helper protocol class that extends the transform routines of the
    Geometry class to a higher-level transformable object, like a Volume
    or Surface. Any class that inherits from Transformable must override
    the geometry() member function for this to work.

    NOTE: For cleanliness, it might be a good idea to actually have Geometry
    inherit from Transformable so that the transform functions are really
    defined here. This is rather cyclical though.
    '''

    def geometry(self):
        '''Returns the geometry associated with the object.'''
        raise NotImplementedError('a class that inherits from Transformable must override geometry()')

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


class Warp:
    '''
    TODOC
    '''

    def __init__(self, warp, source=None, target=None, affine=None):
        self.data = warp
        self.source = source
        self.target = target
        self.affine = affine

    def write(self, filename):
        bindings.morph.write_gcam(self, filename)


def transform_to_midspace(a, b, transform=None, resample=True, interp_method='linear', target=None):
    """
    Transform volumes to midspace.

    Parameters:
        a, b: Images to transform.
        transform: Linear transform aligning volume `a` to `b`.
        resample: Resample images to a single target space. If disabled,
            only the headers are updated.
        interp_method: Interpolation method for resampling.
        target: Resampling target geometry. If None, the target is determined
            by a secondary midspace computed between the *aligned* images. In this
            case, the target shape and resolution is determined by `b`.

    Returns:
        mid_a, mid_b: Midspace-transformed images.
    """
    from scipy.linalg import sqrtm

    if transform is None:
        transform = b.affine @ np.linalg.inv(a.affine)
    else:
        transform = LinearTransform.ensure(transform)
        if transform.type is not None:
            transform = transform.as_ras()
        transform = transform.matrix

    midspace = np.real(sqrtm(transform))
    midspace_inv = np.real(sqrtm(np.linalg.inv(transform)))

    mid_a = a.copy()
    mid_b = b.copy()

    mid_a.affine = midspace @ mid_a.affine
    mid_b.affine = midspace_inv @ mid_b.affine

    if resample:
        if target is None:
            affine = np.real(sqrtm(mid_a.affine @ np.linalg.inv(mid_b.affine))) @ mid_b.affine
            target = Geometry(shape=mid_b.shape, voxsize=mid_b.voxsize, affine=affine)
        mid_a = mid_a.resample_like(target, interp_method=interp_method)
        mid_b = mid_b.resample_like(target, interp_method=interp_method)

    return mid_a, mid_b
