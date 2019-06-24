import numpy as np
import scipy.ndimage


def slicing(start, stop):
    '''
    Creates an N-dimensional slicing from start and stop coordinates.
    '''
    return tuple([slice(int(x), int(y)) for x, y in zip(start, stop)])


def slicing_shape(slices):
    '''
    Array shape defined by a slicing.
    '''
    return tuple([slc.stop - slc.start for slc in slices])


def slicing_origin(slices):
    '''
    Origin (start indices) of a slicing.
    '''
    return tuple([slc.start for slc in slices])


def bbox(mask):
    '''
    Bounding box around an object in a binary image.
    '''
    return scipy.ndimage.find_objects(mask)[0]


def cmass(image):
    '''
    Center of mass of an image.
    '''
    return scipy.ndimage.center_of_mass(image)


def conform(image, shape):
    '''
    Conforms an image to a particular shape. Dimensions will be expanded or cropped
    around the center of the image.
    '''
    diff = (np.array(shape) - image.shape) / 2

    # crop if necessary
    masked = np.abs(diff * (diff < 0))
    offset = np.floor(masked).astype(int)
    conformed = image[slicing(offset, image.shape - np.ceil(masked))]

    # pad if necessary
    masked = diff * (diff > 0)
    padding = np.vstack((np.floor(masked), np.ceil(masked))).T.astype(int)
    conformed = np.pad(conformed, padding, 'constant')

    return conformed, np.floor(masked).astype(int) - offset



def bbox2_3D(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, rmax, cmin, cmax, zmin, zmax


def sample_patch(mri, point, nx, ny, nz, wsize):
    whalf = wsize/2
    patch = np.zeros((wsize,wsize,wsize, mri.shape[3]))
    for xk in range (-whalf,whalf):
        for yk in range (-whalf,whalf):
            for zk in range (-whalf,whalf):
                if (xk == 0 and yk == 0 and zk == 0):
                    xk = 0
                if (xk == 0 and yk == 0 and zk == whalf-1):
                    xk = 0
                xi = int(point[0] + nx[0]*xk + nx[1]*yk + nx[2]*zk)
                yi = int(point[1] + ny[0]*xk + ny[1]*yk + ny[2]*zk)
                zi = int(point[2] + nz[0]*xk + nz[1]*yk + nz[2]*zk)
                val = mri[xi,yi,zi,:]
                patch[xk+whalf,yk+whalf,zk+whalf,:] = val
    return patch
