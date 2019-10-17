import numpy as np
from .. import bindings


class ArrayContainerTemplate:
    '''Internal, abstract template responsible for handling an N-D array. 1D, 2D,
    and 3D implementations represent overlays, images, and volumes, respectively.

    Read and write functions provide default volume-file IO, regardless of the
    dimensionality. This class should only be used as an internal base class, and
    never initialized directly.
    '''

    basedims = None

    def __init__(self, data):
        self.data = np.array(data, copy=False)
        # any extra dimension is assumed to represent data frames
        if self.data.ndim < self.basedims or self.data.ndim > self.basedims + 1:
            raise ValueError('%s (%dD) cannot be initialized by an array with %d dims'
                % (self.__class__.__name__, self.basedims, self.data.ndim))

    @property
    def nframes(self):
        '''Number of data frames.'''
        return self.shape[-1] if self.array.ndim == self.basedims + 1 else 1

    @property
    def shape(self):
        '''Shape of the data array.'''
        return self.data.shape

    @classmethod
    def read(cls, filename):
        result = bindings.vol.read(filename)
        if not isinstance(result, cls):
            warning('reading file "%s" as a %s - not a %s' % (filename, result.__class__.__name__, cls.__name__))
        return result

    @classmethod
    def empty(cls, shape, dtype):
        return cls(np.zeros(shape, dtype, order='F'))


class Overlay(ArrayContainerTemplate):
    '''1D array that represents values corresponding to surface vertices.'''
    basedims = 1


class Image(ArrayContainerTemplate):
    '''2D image with specific geometry.'''
    basedims = 2


class Volume(ArrayContainerTemplate):
    '''3D volume with specific geometry.'''
    basedims = 3

    def __init__(self, data, affine=np.eye(4)):
        super().__init__(data)
        self.affine = affine
        self.voxsize = (1.0, 1.0, 1.0)

    def vox2ras(self):
        return LinearTransform(self.affine)

    def ras2vox(self):
        return self.vox2ras().inverse()

    @property
    def image(self):
        return self.array

    @image.setter
    def image(self, array):
        self.array = array
