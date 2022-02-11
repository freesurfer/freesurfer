import numpy as np
import scipy.ndimage
from collections.abc import Iterable
from scipy.interpolate import interpn

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


def bbox(mask, margin=0):
    '''Bounding box around the object in a binary image.'''
    if not np.any(mask):
        return tuple([slice(0, s) for s in mask.shape])
    bbox = scipy.ndimage.find_objects(mask)[0]
    if margin > 0:
        start = [max(0, c.start - margin) for c in bbox]
        stop = [min(mask.shape[i], c.stop + margin) for i, c in enumerate(bbox)]
        step = [c.step for c in bbox]
        bbox = tuple([slice(*s) for s in zip(start, stop, step)])
    return bbox


def cmass(image):
    '''Center of mass of an image.'''
    return scipy.ndimage.center_of_mass(image)


def resample(source, target_shape, trg2src, interp_method='linear', fill=0, smooth_sigma=0):
    '''
    Resamples a volume array from one space to another given
    a target-to-source transformation matrix.

    Parameters:
        source: Source array to sample from. Must be 3D or 4D.
        target_shape: Shape of the returned target array.
        trg2src: 4x4 affine matrix that transforms target coords to source coords.
        interp_method: Interpolation method. Must be 'linear' or 'nearest'. Default is 'linear'.
        smooth_sigma: Apply gaussian smoothing before resampling (smoothing kernel is
            in voxel space). Default is 0.
    '''
    if source.ndim > 4:
        raise ValueError('resampling can not be done on arrays with more than 4 dimensions')
    elif source.ndim < 3:
        raise NotImplementedError('%dD resampling is not yet supported (must be 3D or 4D)' % source.ndim)
 
    # the resample binding function only works with 4D inputs for easier maintenance, so let's add an axis to any 3D input
    target_shape = target_shape[:3]
    if source.ndim == 4:
        target_shape = (*target_shape, source.shape[-1])

    orig_target_shape = target_shape
    if source.ndim == 3:
        source = source[..., np.newaxis]
        target_shape = (*target_shape, 1)

    if len(target_shape) != source.ndim:
        raise ValueError('resampled target shape (%sD) must match source dims (%sD)' % (len(target_shape), source.ndim))

    if target_shape[-1] != source.shape[-1]:
        raise ValueError('resampled target must have the same number of frames as the source')

    # apply optional gaussian smoothing
    if smooth_sigma > 0:
        sigmas = (smooth_sigma, smooth_sigma, smooth_sigma, 0)
        source = scipy.ndimage.gaussian_filter(source, sigma=sigmas, order=0)  # TODO order?

    trg2src = LinearTransform.ensure(trg2src).matrix

    if interp_method == 'linear':
        return bindings.vol.resample_volume_linear(source, target_shape, trg2src, fill).reshape(orig_target_shape)
    elif interp_method == 'nearest':
        return bindings.vol.resample_volume_nearest(source, target_shape, trg2src, fill).reshape(orig_target_shape)
    else:
        raise ValueError('invalid resample interpolation method "%s"' % interp_method)

def sample_volume(volume, points, interp_method='linear', fill_value=0):
    '''
    Parameters:
        volume: ndarray
            Image array to sample from.
        points: ndarray
            (N, 3) image coordinates to sample from.
        interp_method: str, optional
            The method of interpolation to perform. Supported are 'linear' and 'nearest'.
        fill_value: number, optional
            If provided, the value to use for points outside of the
            interpolation domain. If None, values outside the domain are extrapolated.
    '''
    if volume.ndim > 4:
        raise ValueError('resampling can not be done on arrays with more than 4 dimensions')
    elif volume.ndim < 3:
        raise NotImplementedError('%dD resampling is not yet supported (must be 3D or 4D)' % source.ndim)
    elif volume.ndim == 3:
        volume = volume[..., np.newaxis]

    linvec = [np.arange(0, dim) for dim in volume.shape[:3]]
    sampler = lambda x: interpn(linvec, volume[..., x], points, method=interp_method, bounds_error=False, fill_value=fill_value)
    return np.stack([sampler(f) for f in range(volume.shape[-1])], axis=-1).squeeze()


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


def apply_warp(source, flow, interp_method='linear', indexing='ij'):
    """
    Warps a source image given a flow field. Source and flow volume sizes must match.

    Parameters:
        source: Moving image numpy array. Can be multiframe.
        flow: Flow field numpy array.
        interp_method: Interpolation method. Must be 'linear' or 'nearest'. Default is 'linear'.
        indexing: Cartesian (‘xy’) or matrix (‘ij’) indexing of warp. Default is 'ij'.
    """
    if source.ndim > 4:
        raise ValueError('warping can not be done on arrays with more than 4 dimensions')
    elif source.ndim < 3:
        raise NotImplementedError('%dD warping is not yet supported (must be 3D or 4D)' % source.ndim)

    linvec = [np.arange(0, dim) for dim in source.shape[:3]]
    loc = np.rollaxis(np.array(np.meshgrid(*linvec, indexing=indexing)), 0, 4) + flow
    sample = np.stack((loc[:, :, :, 0], loc[:, :, :, 1], loc[:, :, :, 2]), 3)

    warp = lambda src : interpn(linvec, src, sample, method=interp_method, bounds_error=False, fill_value=0)

    if source.ndim == 3:
        return warp(source)
    if source.ndim == 4:
        return np.stack([warp(source[..., f]) for f in range(source.shape[-1])], axis=-1)
