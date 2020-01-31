import os
import numpy as np
import copy

from . import bindings, warning
from .transform import Transformable, LinearTransform, Geometry


class ArrayContainerTemplate:
    '''
    Internal, abstract template responsible for handling an N-D array. 1D, 2D,
    and 3D implementations represent overlays, images, and volumes, respectively.

    Read and write functions provide default volume-file IO, regardless of the
    dimensionality. This class should only be used as an internal base class, and
    never initialized directly.
    '''

    basedims = None

    def __init__(self, data):
        '''
        Contructs the container object from an array. The input data is not copied, and
        the array should have ndims equal to the subclass' basedims (or basedims + 1).
        Any extra dimension is assumed to represent data frames.
        '''
        # make sure this template class isn't being used directly
        if self.basedims is None:
            raise TypeError('%s should never be initialized directly' % self.__class__.__name)
        self.data = np.array(data, copy=False)
        # extra dim is assumed to represent data frames
        if self.data.ndim < self.basedims or self.data.ndim > self.basedims + 1:
            raise ValueError('%s (%dD) cannot be initialized by an array with %d dims'
                % (self.__class__.__name__, self.basedims, self.data.ndim))

        # any array type might have a valid lookup table
        self.lut = None

    @property
    def nframes(self):
        '''Number of data frames.'''
        return self.shape[-1] if self.data.ndim == self.basedims + 1 else 1

    @property
    def shape(self):
        '''Shape of the data array.'''
        return self.data.shape

    @property
    def dtype(self):
        '''Data type of the array.'''
        return self.data.dtype

    def copy(self):
        '''Returns a deep copy of the instance.'''
        return copy.deepcopy(self)

    @classmethod
    def ensure(cls, unknown):
        '''Ensures that the unknown instance is converted to the correct array class.'''
        if isinstance(unknown, cls):
            return unknown
        if isinstance(unknown, np.ndarray):
            return cls(unknown)
        raise ValueError('Cannot convert object of type %s to %s' % (type(unknown).__name__, cls.__name__))

    @classmethod
    def read(cls, filename):
        '''Reads in array and metadata from volume file.'''
        if not os.path.isfile(filename):
            raise ValueError('file %s does not exist' % filename)
        result = bindings.vol.read(filename)
        # since the volume bindings do all the IO work here, it's possible the returned
        # object type does not match the calling class... if this is the case, print a warning
        if not isinstance(result, cls):
            warning('reading file "%s" as a %s - not a %s' % (filename, result.__class__.__name__, cls.__name__))
        return result

    def write(self, filename):
        '''Writes array and metadata to a volume file.'''
        bindings.vol.write(self, filename)

    @classmethod
    def empty(cls, shape, dtype):
        '''Generates an empty array of given shape and datatype.'''
        return cls(np.zeros(shape, dtype, order='F'))


class Overlay(ArrayContainerTemplate):
    '''1D array that represents values corresponding to surface vertices.'''
    basedims = 1

    def __init__(self, data):
        '''Contructs an overlay from a 1D or 2D data array. The 2nd dimension is
        always assumed to be the number of frames.'''
        super().__init__(data)


class Image(ArrayContainerTemplate, Transformable):
    '''2D image with specific geometry.'''
    basedims = 2

    def __init__(self, data, affine=None):
        '''Contructs an image from a 2D or 3D data array. The 3rd dimension is
        always assumed to be the number of frames.'''
        ArrayContainerTemplate.__init__(self, data)
        self.affine = affine
        self.voxsize = (1.0, 1.0, 1.0)

    def geometry(self):
        '''Returns volume geometry as a `Geometry` instance.'''
        return Geometry(self.shape, self.voxsize, self.affine)


class Volume(ArrayContainerTemplate, Transformable):
    '''
    3D volume with specific geometry.
    
    Attributes:
        data: Pointer to internal 3D (or 4D) array.
        image: Alias to data member (legacy).
        shape: Shape of the internal volume array.
        voxsize: Voxel sizes in millimeters.
        affine: 4x4 vox-to-ras transform matrix.
        nframes: Number of volume frames.
        lut: Embedded label lookup-table for segmentations.
        te: Scan echo time in ms.
        tr: Scan repetition time in ms.
        ti: Scan inversion time in ms.
        flip_angle: Scan flip angle in degrees.
    '''
    basedims = 3

    def __init__(self, data, affine=None, voxsize=None):
        '''
        Contructs a volume from a 3D or 4D data array. The 4th dimension is
        always assumed to be the number of frames.
        '''

        # scan parameters
        self.te = None
        self.tr = None
        self.ti = None
        self.flip_angle = None

        # lookup table for segmentations
        self.lut = None

        # TODEP - this is a temporary fix to support the previous way of loading from a file - it
        # is not an ideal way of handling things and should be removed as soon as possible
        if isinstance(data, str):
            print('moving foward, please load volumes via fs.Volume.read(filename)')
            result = Volume.read(data)
            if not isinstance(result, Volume):
                classname = result.__class__.__name__
                raise ValueError('file "%s" must be read as a %s, like this: %s.read(filename)' % (data, classname, classname))
            data = result.data
            affine = result.affine
            voxsize = result.voxsize
            self.lut = result.lut
            self.te = result.te
            self.tr = result.tr
            self.ti = result.ti
            self.flip_angle = result.flip_angle

        ArrayContainerTemplate.__init__(self, data)
        self.affine = affine
        self.voxsize = voxsize if voxsize else (1.0, 1.0, 1.0)

    def geometry(self):
        '''Returns volume geometry as a `Geometry` instance.'''
        return Geometry(self.shape, self.voxsize, self.affine)

    @property
    def image(self):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        return self.data  # TODEP

    @image.setter
    def image(self, array):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        self.data = array  # TODEP
