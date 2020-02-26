import os
import numpy as np
import copy

from . import bindings, warning
from .geom import resample
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
        # get class name
        classname = self.__class__.__name__

        # make sure this template class isn't being used directly
        if self.basedims is None:
            raise TypeError('%s should never be initialized directly' % classname)

        # make sure a string isn't being provided
        if isinstance(data, str):
            raise ValueError('if loading from file, use the %s.read() class method' % classname)

        self.data = np.array(data, copy=False)
        # extra dim is assumed to represent data frames
        if self.data.ndim < self.basedims or self.data.ndim > self.basedims + 1:
            raise ValueError('%s (%dD) cannot be initialized by an array with %d dims' % (classname, self.basedims, self.data.ndim))

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

        # allow for Volumes to be constructed from arrays with less than 3
        # dimensions in order to act as a catch-all. The appropriate array container
        # classes should really be used if possible since this 'feature' will
        # totally breakdown with multiframe data
        if data.ndim < 3:
            newaxes = [1] * (3 - data.ndim)
            data = data.reshape(*data.shape, *newaxes)

        ArrayContainerTemplate.__init__(self, data)
        self.affine = affine
        self.voxsize = voxsize if voxsize else (1.0, 1.0, 1.0)

    def __getitem__(self, idx):
        '''Returns the cropped volume with a recomputed affine matrix.'''
        return self.crop(idx)

    def geometry(self):
        '''Returns volume geometry as a `Geometry` instance.'''
        return Geometry(self.shape, self.voxsize, self.affine)

    def copy_metadata(self, vol):
        '''Copies metadata from another volume.'''
        self.lut = vol.lut
        self.te = vol.te
        self.tr = vol.tr
        self.ti = vol.ti
        self.flip_angle = vol.flip_angle

    def copy_geometry(self, vol):
        '''Copies voxsize and affine information from another volume.'''
        self.affine = vol.affine
        self.voxsize = vol.voxsize

    def reslice(self, voxsize, smooth_sigma=0):
        '''
        Returns the resampled volume with a given resolution determined by voxel
        size in mm.

        Parameter:
            voxsize: Voxel size of target volume.
            smooth_sigma: Apply gaussian smoothing before resampling (kernel size is
                in voxel space). Default is 0.
        '''
        src_shape = self.shape[:3]
        target_shape = tuple(np.ceil(np.array(self.voxsize).astype(float) * src_shape / voxsize).astype(int))
        
        # get source-to-RAS matrix
        src2ras = self.affine if self.affine is not None else fs.transform.LIA(src_shape, self.voxsize)

        # compute target-to-RAS matrix
        pcrs = np.append(np.array(src_shape) / 2, 1)
        cras = np.matmul(src2ras, pcrs)[:3]
        trg2ras = np.eye(4)
        trg2ras[:3, :3] = self.affine[:3, :3] * voxsize / self.voxsize
        pcrs = np.append(np.array(target_shape) / 2, 1)
        trg2ras[:3, 3] = cras - np.matmul(trg2ras, pcrs)[:3]

        # compute target-to-source matrix
        ras2src = np.linalg.inv(src2ras)
        trg2src = np.matmul(ras2src, trg2ras)

        # resample into new volume
        if self.data.ndim != 3:
            target_shape = (*target_shape, self.nframes)
        resliced_data = resample(self.data, target_shape, trg2src, smooth_sigma=smooth_sigma)
        resliced_vol = Volume(resliced_data, affine=trg2ras, voxsize=voxsize)
        resliced_vol.copy_metadata(self)
        return resliced_vol

    def crop(self, cropping):
        '''
        Returns the cropped volume with a recomputed affine matrix.
        Avoid using this function directly, and instead use direct indexing
        on the volume object, like so:

            cropped = vol[:, 10:-10, :]

        Parameter:
            cropping: Tuple of crop indices (slices).
        '''
        # convert cropping into list
        idxs = [cropping] if not isinstance(cropping, (tuple, list)) else cropping
        
        # we need to extract starting coordinates, and dealing with ellipsis is
        # way too much work for now, so let's not support it
        if Ellipsis in idxs:
            raise NotImplementedError('Ellipsis is not yet supported in cropping indices')
        
        # extract the starting coordinate of the cropping
        start_coords = []
        for i in idxs:
            if isinstance(i, slice):
                start_coords.append(i.start if i.start is not None else 0)
            elif isinstance(i, int):
                raise ValueError('volume cropping does not support dimension removal - crop the data array directly')
            else:
                raise ValueError('volume cropping must be indexed by slices (:) only')

        # crop the raw array
        cropped_data = self.data[cropping]

        # compute new affine if one exists
        if self.affine is not None:
            matrix = np.eye(4)
            matrix[:3, :3] = self.affine[:3, :3]
            p0 = self.vox2ras().transform(start_coords[:3])
            matrix[:3, 3] = p0
            pcrs = np.append(np.array(cropped_data.shape[:3]) / 2, 1)
            cras = np.matmul(matrix, pcrs)[:3]
            matrix[:3, 3] = 0
            matrix[:3, 3] = cras - np.matmul(matrix, pcrs)[:3]
        else:
            matrix = None

        # construct cropped volume
        cropped_vol = Volume(cropped_data, affine=matrix, voxsize=self.voxsize)
        cropped_vol.copy_metadata(self)
        return cropped_vol

    @property
    def image(self):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        return self.data  # TODEP

    @image.setter
    def image(self, array):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        self.data = array  # TODEP
