import numpy as np
import scipy.ndimage
from collections.abc import Iterable

from . import bindings
from .transform import LinearTransform


class Slicing(tuple):
    '''Slice tuple used for indexing subregions of a numpy array.'''

    def __new__(cls, start, stop):
        if len(start) != len(stop):
            raise ValueError('start coord (%dD) does not match stop coord (%dD))' % (len(start), len(stop)))
        return super(Slicing, cls).__new__(cls, [slice(s, t) for s, t in zip(start, stop)])

    @property
    def start(self):
        return tuple([s.start for s in self])

    @property
    def stop(self):
        return tuple([s.stop for s in self])

    @property
    def shape(self):
        return tuple([s.stop - s.start for s in self])

    def grow(self, dist):
        if not isinstance(dist, Iterable):
            dist = [dist] * len(self)
        start = [s.start - int(d) for s, d in zip(self, dist)]
        stop  = [s.stop  + int(d) for s, d in zip(self, dist)]
        return Slicing(start, stop)

    def shrink(self, dist):
        if isinstance(dist, Iterable):
            dist = np.array(dist)
        return self.grow(dist * -1)


def bbox(mask):
    '''Bounding box around the object in a binary image.'''
    return scipy.ndimage.find_objects(mask)[0]


def cmass(image):
    '''Center of mass of an image.'''
    return scipy.ndimage.center_of_mass(image)


def resample(source, target_shape, trg2src, smooth_sigma=0):
    '''
    Resamples a volume array from one space to another given
    a target-to-source transformation matrix.

    Parameters:
        source: Source array to sample from. Must be 3D or 4D.
        target_shape: Shape of the returned target array.
        trg2src: 4x4 affine matrix that transforms target coords to source coords.
        smooth_sigma: Apply gaussian smoothing before resampling. Default is 0.
    '''
    if source.ndim > 4:
        raise ValueError('resampling can not be done on arrays with more than 4 dimensions')
    elif source.ndim < 3:
        raise NotImplementedError('%dD resampling is not yet supported (must be 3D or 4D)' % source.ndim)

    if len(target_shape) != source.ndim:
        raise ValueError('resampled target shape (%sD) must match source dims (%sD)' % (len(target_shape), source.ndim))
 
    # the resample binding function only works with 4D inputs for easier maintenance, so let's add an axis to any 3D input
    orig_target_shape = target_shape
    if source.ndim == 3:
        source = source[..., np.newaxis]
        target_shape = (*target_shape, 1)

    if target_shape[-1] != source.shape[-1]:
        raise ValueError('resampled target must have the same number of frames as the source')

    # apply optional gaussian smoothing
    if smooth_sigma > 0:
        sigmas = (smooth_sigma, smooth_sigma, smooth_sigma, 0)
        source = scipy.ndimage.gaussian_filter(source, sigma=sigmas, order=0)  # TODO order?

    trg2src = LinearTransform.ensure(trg2src).matrix
    return bindings.vol.resample_volume(source, target_shape, trg2src).reshape(orig_target_shape)


def sample_into_volume(volume, weights, coords, values):
    '''
      Interpolates and adds values at given coordinates into a volume and corresponding
      weights array. Since this function performs in-place modification of the volume and
      weights arrays, they MUST be double or float precision to avoid being copied.

      Parameters:
        volume: Float or double array to add sampled values into.
        weights: Float or double array to add sampling weights into.
        coords: List of volume sampling coordinates. 
        values: List of values at each coordinate.
    '''
    if volume.ndim > 4:
        raise ValueError('resampling can not be done on arrays with more than 4 dimensions')
    elif volume.ndim < 3:
        raise NotImplementedError('%dD resampling is not yet supported (must be 3D or 4D)' % source.ndim)
    bindings.vol.sample_into_volume(volume, weights, coords, values)
