import os
import numpy as np
import copy
import scipy

from collections.abc import Iterable

from . import bindings, warning
from .geom import resample, apply_warp
from .transform import Transformable, LinearTransform, Geometry, LIA
from . import orientation as otn


class ArrayContainerTemplate:
    '''
    Internal, abstract template responsible for handling an N-D array. 1D, 2D,
    and 3D implementations represent overlays, images, and volumes, respectively.

    Read and write functions provide default volume-file IO, regardless of the
    dimensionality. This class should only be used as an internal base class, and
    never initialized directly.
    '''

    basedims = None

    def __init__(self, data, lut=None, swap_batch_dim=False):
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

        # swap frame dimension from first axis to last
        if swap_batch_dim and data.ndim > self.basedims:
            data = np.moveaxis(data, 0, -1)

        self.data = np.array(data, copy=False)
        # extra dim is assumed to represent data frames
        if self.data.ndim < self.basedims or self.data.ndim > self.basedims + 1:
            raise ValueError('%s (%dD) cannot be initialized by an array with %d dims' % (classname, self.basedims, self.data.ndim))

        # any array type might have a valid lookup table
        self.lut = lut

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

    def copy(self, data=None):
        '''
        Returns a deep copy of the instance.

        Parameters:
            data: Replace internal data array. Default is None.
        '''
        if data is not None:
            # to save memory if we're replacing the data, make sure not to create
            # a duplicate instance of the original data before replacing it
            shallow = copy.copy(self)
            shallow.data = None
            copied = copy.deepcopy(shallow)
            copied.data = data
            return copied
        else:
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

    def __init__(self, data, lut=None, **kwargs):
        '''Contructs an overlay from a 1D or 2D data array. The 2nd dimension is
        always assumed to be the number of frames.'''
        super().__init__(data, lut=lut, **kwargs)


class Image(ArrayContainerTemplate, Transformable):
    '''2D image with specific geometry.'''
    basedims = 2

    def __init__(self, data, affine=None, pixsize=None, lut=None, **kwargs):
        '''Contructs an image from a 2D or 3D data array. The 3rd dimension is
        always assumed to be the number of frames.'''
        ArrayContainerTemplate.__init__(self, data, lut=lut, **kwargs)
        self.affine = affine
        self.pixsize = pixsize if pixsize is not None else (1.0, 1.0)

    def geometry(self):
        '''Returns volume geometry as a `Geometry` instance.'''
        return Geometry(self.shape, self.pixsize, self.affine)

    def copy_geometry(self, image):
        '''Copies pixsize and affine information from another image.'''
        self.affine = image.affine
        self.pixsize = image.pixsize

    def reslice(self, pixsize, interp_method='linear', smooth_sigma=0):
        '''
        Returns the resampled image with a given resolution determined by pixel size in mm.

        Parameters:
            pixsize: Voxel size of target volume. Can be single value or list.
            interp_method: Interpolation method. Must be 'linear' or 'nearest'. Default is 'linear'.
            smooth_sigma: Apply gaussian smoothing before resampling (kernel size is
                in voxel space). Default is 0.
        '''

        # convert single value to list
        if not isinstance(pixsize, Iterable):
            pixsize = [pixsize] * 2

        # reshape image into '3D' array
        src_shape_3d = (*self.shape[:2], 1)
        if self.data.ndim != 2:
            src_shape_3d = (*src_shape_3d, self.nframes)

        # reslice the image as a 'volume'
        src_vol = Volume(np.reshape(self.data, src_shape_3d), affine=self.affine, voxsize=(*self.pixsize[:2], 1.0))
        trg_vol = src_vol.reslice((*pixsize[:2], 1.0), interp_method=interp_method, smooth_sigma=smooth_sigma)

        # transfer back into image
        trg_shape = trg_vol.shape[:2]
        if self.data.ndim != 2:
            trg_shape = (*trg_shape, self.nframes)
        resliced_image = Image(np.reshape(trg_vol.data, trg_shape), affine=trg_vol.affine, pixsize=trg_vol.voxsize[:2])
        return resliced_image


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

    def __init__(self, data, affine=None, voxsize=None, lut=None, **kwargs):
        '''
        Contructs a volume from a 3D or 4D data array. The 4th dimension is
        always assumed to be the number of frames.
        '''

        # scan parameters
        self.te = None
        self.tr = None
        self.ti = None
        self.flip_angle = None

        # allow for Volumes to be constructed from arrays with less than 3
        # dimensions in order to act as a catch-all. The appropriate array container
        # classes should really be used if possible since this 'feature' will
        # totally breakdown with multiframe data
        if data.ndim < 3:
            newaxes = [1] * (3 - data.ndim)
            data = data.reshape(*data.shape, *newaxes)

        ArrayContainerTemplate.__init__(self, data, lut=lut, **kwargs)
        self.affine = affine
        self.voxsize = voxsize if voxsize is not None else (1.0, 1.0, 1.0)

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

    def reslice(self, voxsize, interp_method='linear', smooth_sigma=0):
        '''
        Returns the resampled volume with a given resolution determined by voxel
        size in mm.

        Parameters:
            voxsize: Voxel size of target volume. Can be single value or list.
            interp_method: Interpolation method. Must be 'linear' or 'nearest'. Default is 'linear'.
            smooth_sigma: Apply gaussian smoothing before resampling (kernel size is
                in voxel space). Default is 0.
        '''

        # convert single value to list
        if not isinstance(voxsize, Iterable):
            voxsize = [voxsize] * self.basedims

        if np.allclose(voxsize, self.voxsize):
            return self

        src_shape = self.shape[:3]
        target_shape = tuple(np.ceil(np.array(self.voxsize).astype(float) * src_shape / voxsize).astype(int))
        
        # get source-to-RAS matrix
        src2ras = self.affine if self.affine is not None else LIA(src_shape, self.voxsize)

        # compute target-to-RAS matrix
        pcrs = np.append(np.array(src_shape) / 2, 1)
        cras = np.matmul(src2ras, pcrs)[:3]
        trg2ras = np.eye(4)
        trg2ras[:3, :3] = src2ras[:3, :3] * voxsize / self.voxsize
        pcrs = np.append(np.array(target_shape) / 2, 1)
        trg2ras[:3, 3] = cras - np.matmul(trg2ras, pcrs)[:3]

        # compute target-to-source matrix
        ras2src = np.linalg.inv(src2ras)
        trg2src = np.matmul(ras2src, trg2ras)

        # resample into new volume
        if self.data.ndim != 3:
            target_shape = (*target_shape, self.nframes)
        resliced_data = resample(self.data, target_shape, trg2src, interp_method=interp_method, smooth_sigma=smooth_sigma)
        resliced_vol = Volume(resliced_data, affine=trg2ras, voxsize=voxsize)
        resliced_vol.copy_metadata(self)
        return resliced_vol

    def crop(self, cropping):
        '''
        Returns the cropped volume with a recomputed affine matrix.
        Avoid using this function directly, and instead use direct indexing
        on the volume object, like so:

            cropped = vol[:, 10:-10, :]

        Parameters:
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

    def bbox(self, thresh=0, margin=0):
        '''
        TODOC
        '''
        mask = self.data > thresh
        if not np.any(mask):
            return tuple([slice(0, s) for s in mask.shape])
        cropping = scipy.ndimage.find_objects(mask)[0]
        if margin > 0:
            start = [max(0, c.start - margin) for c in cropping]
            stop = [min(self.shape[i], c.stop + margin) for i, c in enumerate(cropping)]
            step = [c.step for c in cropping]
            cropping = tuple([slice(*s) for s in zip(start, stop, step)])
        return cropping

    def crop_to_bbox(self, thresh=0, margin=0):
        '''
        TODOC
        '''
        cropping = self.bbox(thresh=thresh, margin=margin)
        return self[cropping]

    def fit_to_shape(self, shape, center='image'):
        '''
        Returns a volume fit to a given shape. Image will be
        centered in the conformed volume.

        TODO: Enable multi-frame support.
        '''

        # This is a quick hack
        if center == 'bbox':
            v = self.crop_to_bbox()
            return v.fit_to_shape(shape, center='image')

        if self.nframes > 1:
            raise NotImplementedError('multiframe volumes not support yet for shape refit')

        delta = (np.array(shape) - np.array(self.shape[:3])) / 2
        low = np.floor(delta).astype(int)
        high = np.ceil(delta).astype(int)

        c_low = np.clip(low, 0, None)
        c_high = np.clip(high, 0, None)
        conformed_data = np.pad(self.data.squeeze(), list(zip(c_low, c_high)), mode='constant')

        # note: low and high are intentionally swapped here
        c_low = np.clip(-high, 0, None)
        c_high = conformed_data.shape[:3] - np.clip(-low, 0, None)
        cropping = tuple([slice(a, b) for a, b in zip(c_low, c_high)])
        conformed_data = conformed_data[cropping]

        # compute new affine if one exists
        if self.affine is not None:
            matrix = np.eye(4)
            matrix[:3, :3] = self.affine[:3, :3]
            p0crs = np.clip(-high, 0, None) - np.clip(low, 0, None)
            p0 = self.vox2ras().transform(p0crs)
            matrix[:3, 3] = p0
            pcrs = np.append(np.array(conformed_data.shape[:3]) / 2, 1)
            cras = np.matmul(matrix, pcrs)[:3]
            matrix[:3, 3] = 0
            matrix[:3, 3] = cras - np.matmul(matrix, pcrs)[:3]
        else:
            matrix = None

        # construct cropped volume
        conformed_vol = Volume(conformed_data, affine=matrix, voxsize=self.voxsize)
        conformed_vol.copy_metadata(self)
        return conformed_vol

    def transform(self, trf, interp_method='linear', indexing='ij'):
        '''
        Returns a volume transformed by either a dense deformation field or affine
        target-to-source matrix (4x4).

        Parameters:
            trf: Transform to apply. Can be a dense warp or an affine matrix.
            interp_method: Interpolation method. Must be 'linear' or 'nearest'. Default is 'linear'.
            indexing: Cartesian (‘xy’) or matrix (‘ij’) indexing of warp. Default is 'ij'.
        '''

        # convert high-level types to numpy arrays
        trf_target = None
        if isinstance(trf, Volume):
            trf = trf.data
        elif isinstance(trf, LinearTransform):
            trf_target = trf.target
            trf = trf.matrix

        if trf_target is None:
            trf_target = self

        # assert transform type and apply
        if trf.shape[-1] == self.basedims:
            resampled_data = apply_warp(self.data, trf, interp_method=interp_method, indexing=indexing)
        elif len(trf.shape) == 2:
            resampled_data = resample(self.data, self.shape, trf, interp_method=interp_method)
        else:
            raise ValueError('cannot determine transform type')

        # construct new volume
        resampled = Volume(resampled_data)
        resampled.copy_geometry(trf_target)
        resampled.copy_metadata(self)
        return resampled

    def reorient(self, orientation):
        """
        Realigns image data and world matrix to conform to a specific slice orientation.

        TODO ensure header if coorectly updated for mutlires data
        """
        trg_orientation = orientation.upper()
        src_orientation = otn.orientation_from_matrix(self.affine)
        if trg_orientation == src_orientation.upper():
            return self

        # extract world axes
        get_world_axes = lambda aff: np.argmax(np.absolute(np.linalg.inv(aff)), axis=0)
        trg_matrix = otn.matrix_from_orientation(trg_orientation)
        src_matrix = otn.matrix_from_orientation(src_orientation)
        world_axes_trg = get_world_axes(trg_matrix[:self.basedims, :self.basedims])
        world_axes_src = get_world_axes(src_matrix[:self.basedims, :self.basedims])

        voxsize = np.asarray(self.voxsize)
        voxsize = voxsize[world_axes_src][world_axes_trg]

        # init
        data = self.data.copy()
        affine = self.affine.copy()

        # align axes
        affine[:, world_axes_trg] = affine[:, world_axes_src]
        for i in range(self.basedims):
            if world_axes_src[i] != world_axes_trg[i]:
                data = np.swapaxes(data, world_axes_src[i], world_axes_trg[i])
                swapped_axis_idx = np.where(world_axes_src == world_axes_trg[i])
                world_axes_src[swapped_axis_idx], world_axes_src[i] = world_axes_src[i], world_axes_src[swapped_axis_idx]

        # align directions
        dot_products = np.sum(affine[:3, :3] * trg_matrix[:3, :3], axis=0)
        for i in range(self.basedims):
            if dot_products[i] < 0:
                data = np.flip(data, axis=i)
                affine[:, i] = - affine[:, i]
                affine[:3, 3] = affine[:3, 3] - affine[:3, i] * (data.shape[i] - 1)

        reoriented = Volume(data, affine, voxsize=voxsize)
        reoriented.copy_metadata(self)
        return reoriented

    def conform(self, shape=None, voxsize=1.0, orientation='LIA', interp_method='linear', 
                dtype=None, smooth_sigma=0):
        """
        Conforms image to a specific shape, type, resolution, and orientation.
        """
        conformed = self.reorient(orientation)
        conformed = conformed.reslice(voxsize, interp_method=interp_method, 
                                      smooth_sigma=smooth_sigma)
        if shape is not None:
            conformed = conformed.fit_to_shape(shape)
        if dtype is not None:
            conformed.data = conformed.data.astype(dtype)
        return conformed

    def resample_like(self, target, interp_method='linear'):
        '''
        Returns a resampled image in the target space.
        '''
        if target.affine is None or self.affine is None:
            raise ValueError("Can't resample volume without geometry information.")

        vox2vox = LinearTransform.matmul(self.ras2vox(), target.vox2ras())
        resampled_data = resample(self.data, target.shape, vox2vox, interp_method=interp_method)

        resampled = Volume(resampled_data)
        resampled.copy_geometry(target)
        resampled.copy_metadata(self)
        return resampled

    @property
    def image(self):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        return self.data  # TODEP

    @image.setter
    def image(self, array):
        '''Internal data array. This will be deprecated - just used the Volume.data member.'''
        self.data = array  # TODEP
