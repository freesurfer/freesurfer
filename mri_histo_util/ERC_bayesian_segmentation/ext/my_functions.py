import os
import numpy as np
import nibabel as nib
import cv2
import torch
import scipy.ndimage
import scipy.sparse as sp
import ext.interpol as interpol
from torch.nn import functional
from torch.utils.data import Dataset, DataLoader

###############################3

# Crop label volume
def cropLabelVol(V,
                 margin=10,
                 threshold=0):

    # Make sure it's 3D
    margin = np.array(margin)
    if len(margin.shape) < 2:
        margin = [margin, margin, margin]

    if len(V.shape) < 2:
        V = V[..., np.newaxis]
    if len(V.shape) < 3:
        V = V[..., np.newaxis]

    # Now
    idx = np.where(V > threshold)
    i1 = np.max([0, np.min(idx[0]) - margin[0]]).astype('int')
    j1 = np.max([0, np.min(idx[1]) - margin[1]]).astype('int')
    k1 = np.max([0, np.min(idx[2]) - margin[2]]).astype('int')
    i2 = np.min([V.shape[0], np.max(idx[0]) + margin[0] + 1]).astype('int')
    j2 = np.min([V.shape[1], np.max(idx[1]) + margin[1] + 1]).astype('int')
    k2 = np.min([V.shape[2], np.max(idx[2]) + margin[2] + 1]).astype('int')

    cropping = [i1, j1, k1, i2, j2, k2]
    cropped = V[i1:i2, j1:j2, k1:k2]

    return cropped, cropping

###############################3

def applyCropping(V, cropping):
    i1 = cropping[0]
    j1 = cropping[1]
    k1 = cropping[2]
    i2 = cropping[3]
    j2 = cropping[4]
    k2 = cropping[5]

    if len(V.shape)>2:
        Vcropped = V[i1:i2, j1: j2, k1: k2, ...]
    else:
        Vcropped = V[i1:i2, j1: j2]

    return Vcropped

###############################3

def viewVolume(x, aff=None):

    if aff is None:
        aff = np.eye(4)
    else:
        if type(aff) == torch.Tensor:
            aff = aff.detach().cpu().numpy()

    if type(x) is not list:
        x = [x]

    cmd = 'source /usr/local/freesurfer/nmr-dev-env-bash && freeview '

    for n in np.arange(len(x)):
        vol = x[n]
        if type(vol) == torch.Tensor:
            vol = vol.detach().cpu().numpy()
        vol = np.squeeze(np.array(vol))
        name = '/tmp/' + str(n) + '.nii.gz'
        MRIwrite(vol, aff, name)
        cmd = cmd + ' ' + name

    os.system(cmd + ' &')

###############################3

def MRIwrite(volume, aff, filename, dtype=None):

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)

    nib.save(nifty, filename)

###############################3

def MRIread(filename, dtype=None, im_only=False):

    assert filename.endswith(('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    volume = x.get_fdata()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

###############################3

def downsampleMRI2d(X, aff, shape, factors, mode='image'):

    assert False, 'Function not debugged/tested yet...'

    assert mode=='image' or mode=='labels', 'Mode must be image or labels'
    assert (shape is None) or (factors is None), 'Either shape or factors must be None'
    assert (shape is not None) or (factors is not None), 'Either shape or factors must be not None'

    if shape is not None:
        factors = np.array(shape) / X.shape[0:2]
    else:
        factors = np.array(factors)
        shape = np.round(X.shape[0:2] * factors).astype('int')

    if mode == 'image':
        if np.mean(factors) < 1: # shrink
            Y = cv2.resize(X, shape, interpolation=cv2.INTER_AREA)
        else:  # expan
            Y = cv2.resize(X, shape, interpolation=cv2.INTER_LINEAR)
    else:
        Y = cv2.resize(X, shape, interpolation=cv2.INTER_NEAREST)

    aff2 = aff
    aff2[:, 0] = aff2[:, 0] * factors[0]
    aff2[:, 1] = aff2[:, 1] * factors[1]
    aff2[0:3, 3] = aff2[0:3, 3] + aff[0:3, 0:3] * (0.5*np.array([[factors[0]], [factors[1]], [1]])-0.5)

    return Y, aff2

###############################3

def vox2ras(vox, vox2ras):

    vox2 = np.concatenate([vox, np.ones(shape=[1, vox.shape[1]])], axis=0)

    ras = np.matmul(vox2ras, vox2)[:-1, :]

    return ras

###############################

def ras2vox(ras, vox2ras):

    ras2 = np.concatenate([ras, np.ones(shape=[1, ras.shape[1]])], axis=0)

    vox = np.matmul(np.linalg.inv(vox2ras), ras2)[:-1, :]

    return vox


###############################3

def prepBiasFieldBase3d(siz, max_order):
    x = np.linspace(-1, 1, siz[0])
    y = np.linspace(-1, 1, siz[1])
    z = np.linspace(-1, 1, siz[2])
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    PSI = []
    for o in range(max_order + 1):
        for ox in range(o + 1):
            for oy in range(o + 1):
                for oz in range(o + 1):
                    if (ox + oy + oz) == o:
                        psi = np.ones(siz)
                        for i in range(1, ox + 1):
                            psi = psi * xx
                        for j in range(1, oy + 1):
                            psi = psi * yy
                        for k in range(1, oz + 1):
                            psi = psi * zz
                        PSI.append(psi)

    PSI = np.stack(PSI, axis=-1)

    return PSI

###############################3

def grad3d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3, 1]))
    Gz = scipy.ndimage.convolve(X, np.reshape(h, [1, 1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy + Gz * Gz)

    if provide_gradients:
        return Gmodule, Gx, Gy, Gz
    else:
        return Gmodule

###############################3

def grad2d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy)

    if provide_gradients:
        return Gmodule, Gx, Gy
    else:
        return Gmodule


########################
def torch_resize(I, aff, resolution, device, power_factor_at_half_width=5, dtype=torch.float32, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_resize(I, aff, resolution, device, power_factor_at_half_width, dtype, slow)

    slow = slow or device.type == 'cpu'
    voxsize = np.sqrt(np.sum(aff[:-1, :-1] ** 2, axis=0))
    newsize = np.round(I.shape[0:3] * (voxsize / resolution)).astype(int)
    factors = np.array(I.shape[0:3]) / np.array(newsize)
    k = np.log(power_factor_at_half_width) / np.pi
    sigmas = k * factors
    sigmas[sigmas<=k] = 0  # TODO: we could maybe remove this line, to make sure we always smooth a bit?

    if len(I.shape) not in (3, 4):
        raise Exception('torch_resize works with 3D or 3D+label volumes')
    no_channels = len(I.shape) == 3
    if no_channels:
        I = I[:, :, :, None]
    if torch.is_tensor(I):
        I = I.permute([3, 0, 1, 2])
    else:
        I = I.transpose([3, 0, 1, 2])

    It_lowres = None
    for c in range(len(I)):
        It = torch.as_tensor(I[c], device=device, dtype=dtype)[None, None]
        # Smoothen if needed
        for d in range(3):
            if sigmas[d]>0:
                sl = np.ceil(sigmas[d] * 2.5).astype(int)
                v = np.arange(-sl, sl + 1)
                gauss = np.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / np.sum(gauss)
                kernel = torch.tensor(kernel,  device=device, dtype=dtype)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None,  None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])

            It = It.permute([0, 1, 4, 2, 3])
        It = torch.squeeze(It)
        It, aff2 = myzoom_torch(It, aff, newsize, device)
        It = It.detach()
        if torch.is_tensor(I):
            It = It.to(I.device)
        else:
            It = It.cpu().numpy()
        if len(I) == 1:
            It_lowres = It[None]
        else:
            if It_lowres is None:
                if torch.is_tensor(It):
                    It_lowres = It.new_empty([len(I), *It.shape])
                else:
                    It_lowres = np.empty_like(It, shape=[len(I), *It.shape])
            It_lowres[c] = It

        torch.cuda.empty_cache()

    if not no_channels:
        if torch.is_tensor(I):
            It_lowres = It_lowres.permute([1, 2, 3, 0])
        else:
            It_lowres = It_lowres.transpose([1, 2, 3, 0])
    else:
        It_lowres = It_lowres[0]

    return It_lowres, aff2


@torch.jit.script
def conv_slow_fallback(x, kernel):
    """1D Conv along the last dimension with padding"""
    y = torch.zeros_like(x)
    x = torch.nn.functional.pad(x, [(len(kernel) - 1) // 2]*2)
    x = x.unfold(-1, size=len(kernel), step=1)
    x = x.movedim(-1, 0)
    for i in range(len(kernel)):
        y = y.addcmul_(x[i], kernel[i])
    return y


########################

def myzoom_torch(X, aff, newsize, device):

    if len(X.shape)==3:
        X = X[..., None]

    factors = np.array(newsize) / np.array(X.shape[:-1])
    delta = (1.0 - factors) / (2.0 * factors)

    vx = torch.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0], dtype=torch.float, device=device)[:newsize[0]]
    vy = torch.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1], dtype=torch.float, device=device)[:newsize[1]]
    vz = torch.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2], dtype=torch.float, device=device)[:newsize[2]]

    vx[vx < 0] = 0
    vy[vy < 0] = 0
    vz[vz < 0] = 0
    vx[vx > (X.shape[0]-1)] = (X.shape[0]-1)
    vy[vy > (X.shape[1] - 1)] = (X.shape[1] - 1)
    vz[vz > (X.shape[2] - 1)] = (X.shape[2] - 1)

    fx = torch.floor(vx).int()
    cx = fx + 1
    cx[cx > (X.shape[0]-1)] = (X.shape[0]-1)
    wcx = vx - fx
    wfx = 1 - wcx

    fy = torch.floor(vy).int()
    cy = fy + 1
    cy[cy > (X.shape[1]-1)] = (X.shape[1]-1)
    wcy = vy - fy
    wfy = 1 - wcy

    fz = torch.floor(vz).int()
    cz = fz + 1
    cz[cz > (X.shape[2]-1)] = (X.shape[2]-1)
    wcz = vz - fz
    wfz = 1 - wcz

    Y = torch.zeros([newsize[0], newsize[1], newsize[2], X.shape[3]], dtype=torch.float, device=device)

    dtype = X.dtype
    for channel in range(X.shape[3]):
        Xc = X[:,:,:,channel]

        tmp1 = torch.zeros([newsize[0], Xc.shape[1], Xc.shape[2]], dtype=dtype, device=device)
        for i in range(newsize[0]):
            tmp1[i, :, :] = wfx[i] * Xc[fx[i], :, :] +  wcx[i] * Xc[cx[i], :, :]
        tmp2 = torch.zeros([newsize[0], newsize[1], Xc.shape[2]], dtype=dtype, device=device)
        for j in range(newsize[1]):
            tmp2[:, j, :] = wfy[j] * tmp1[:, fy[j], :] +  wcy[j] * tmp1[:, cy[j], :]
        for k in range(newsize[2]):
            Y[:, :, k, channel] = wfz[k] * tmp2[:, :, fz[k]] +  wcz[k] * tmp2[:, :, cz[k]]

    if Y.shape[3] == 1:
        Y = Y[:,:,:, 0]

    if aff is not None:
        aff_new = aff.copy()
        for c in range(3):
            aff_new[:-1, c] = aff_new[:-1, c] / factors[c]
        aff_new[:-1, -1] = aff_new[:-1, -1] - aff[:-1, :-1] @ (0.5 - 0.5 / factors)
        return Y, aff_new
    else:
        return Y


def read_LUT(LUT_file):

    names = [None] * 10000
    colors = [None] * 10000
    with open(LUT_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        split = line.lstrip().split()
        if split and not split[0].startswith('#'):
            index, name = split[:2]
            color = np.asarray(list(map(int, split[2:5])), dtype=np.float64)
            idx = int(index)
            names[idx] = name
            colors[idx] = color

    return names, colors

####

def flip_left_right_labels(X):

    lut = np.zeros(3000)

    lut[2] = 41; lut[4] = 43; lut[5] = 44; lut[7] = 46; lut[8] = 47; lut[10] = 49; lut[11] = 50; lut[12] = 51; lut[13] = 52; lut[14] = 14;
    lut[15] = 15; lut[16] = 16; lut[17] = 53; lut[18] = 54; lut[24] = 24; lut[26] = 58; lut[28] = 60; lut[41] = 2; lut[43] = 4; lut[44] = 5;
    lut[46] = 7; lut[47] = 8; lut[49] = 10; lut[50] = 11; lut[51] = 12; lut[52] = 13; lut[53] = 17; lut[54] = 18; lut[58] = 26; lut[60] = 28;
    lut[1000:2000] = range(2000, 3000)
    lut[2000:3000] = range(1000, 2000)

    Y = lut[X.astype(int)]

    return Y

def cumprod(sequence, reverse=False, exclusive=False):
    """Perform the cumulative product of a sequence of elements.

    Parameters
    ----------
    sequence : any object that implements `__iter__`
        Sequence of elements for which the `__mul__` operator is defined.
    reverse : bool, default=False
        Compute cumulative product from right-to-left:
        `cumprod([a, b, c], reverse=True) -> [a*b*c, b*c, c]`
    exclusive : bool, default=False
        Exclude self from the cumulative product:
        `cumprod([a, b, c], exclusive=True) -> [1, a, a*b]`

    Returns
    -------
    product : list
        Product of the elements in the sequence.

    """
    if reverse:
        sequence = reversed(sequence)
    accumulate = None
    seq = [1] if exclusive else []
    for elem in sequence:
        if accumulate is None:
            accumulate = elem
        else:
            accumulate = accumulate * elem
        seq.append(accumulate)
    if exclusive:
        seq = seq[:-1]
    if reverse:
        seq = list(reversed(seq))
    return seq


def sub2ind(subs, shape, out=None):
    """Convert sub indices (i, j, k) into linear indices.

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    subs : (D, ...) tensor
        List of sub-indices. The first dimension is the number of dimension.
        Each element should have the same number of elements and shape.
    shape : (D,) vector_like
        Size of each dimension. Its length should be the same as the
        first dimension of ``subs``.
    out : tensor, optional
        Output placeholder

    Returns
    -------
    ind : (...) tensor
        Linear indices
    """
    *subs, ind = subs
    if out is None:
        ind = ind.clone()
    else:
        out.reshape(ind.shape).copy_(ind)
        ind = out
    backend = dict(dtype=ind.dtype, device=ind.device)
    stride = cumprod(shape[1:], reverse=True)
    for i, s in zip(subs, stride):
        ind += torch.as_tensor(i, **backend) * torch.as_tensor(s, **backend)
    return ind


def ind2sub(ind, shape, out=None):
    """Convert linear indices into sub indices (i, j, k).

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    ind : tensor
        Linear indices
    shape : (D,) vector_like
        Size of each dimension.
    out : tensor, optional
        Output placeholder

    Returns
    -------
    subs : (D, ...) tensor
        Sub-indices.
    """
    backend = dict(dtype=ind.dtype, device=ind.device)
    stride = cumprod(shape, reverse=True, exclusive=True)
    stride = torch.as_tensor(stride, **backend)
    if out is None:
        sub = ind.new_empty([len(shape), *ind.shape])
    else:
        sub = out.reshape([len(shape), *ind.shape])
    sub[:, ...] = ind
    for d in range(len(shape)):
        if d > 0:
            torch.remainder(sub[d], torch.as_tensor(stride[d-1], **backend), out=sub[d])
        sub[d] = torch.div(sub[d], stride[d], out=sub[d], rounding_mode='trunc')
    return sub


def affine_sub(affine, shape, indices):
    """Update an affine matrix according to a sub-indexing of the lattice.

    Notes
    -----
    .. Only sub-indexing that *keep an homogeneous voxel size* are allowed.
       Therefore, indices must be `None` or of type `int`, `slice`, `ellipsis`.

    Parameters
    ----------
    affine : (..., ndim_out[+1], ndim_in+1) tensor
        Input affine matrix.
    shape : (ndim_in,) sequence[int]
        Input shape.
    indices : tuple[slice or ellipsis]
        Subscripting indices.

    Returns
    -------
    affine : (..., ndim_out[+1], ndim_new+1) tensor
        Updated affine matrix.
    shape : (ndim_new,) tuple[int]
        Updated shape.

    """
    def is_int(elem):
        if torch.is_tensor(elem):
            return elem.dtype in (torch.int32, torch.int64)
        elif isinstance(elem, int):
            return True
        else:
            return False

    def to_int(elem):
        if torch.is_tensor(elem):
            return elem.item()
        else:
            assert isinstance(elem, int)
            return elem

    # check types
    nb_dim = affine.shape[-1] - 1
    backend = dict(dtype=affine.dtype, device=affine.device)
    if torch.is_tensor(shape):
        shape = shape.tolist()
    if len(shape) != nb_dim:
        raise ValueError('Expected shape of length {}. Got {}'
                         .format(nb_dim, len(shape)))
    if not isinstance(indices, tuple):
        raise TypeError('Indices should be a tuple.')
    indices = list(indices)

    # compute the number of input dimension that correspond to each index
    #   > slice index one dimension but eliipses index multiple dimension
    #     and their number must be computed.
    nb_dims_in = []
    ind_ellipsis = None
    for n_ind, ind in enumerate(indices):
        if isinstance(ind, slice):
            nb_dims_in.append(1)
        elif ind is Ellipsis:
            if ind_ellipsis is not None:
                raise ValueError('Cannot have more than one ellipsis.')
            ind_ellipsis = n_ind
            nb_dims_in.append(-1)
        elif is_int(ind):
            nb_dims_in.append(1)
        elif ind is None:
            nb_dims_in.append(0)
        else:
            raise TypeError('Indices should be None, integers, slices or '
                            'ellipses. Got {}.'.format(type(ind)))
    nb_known_dims = sum(nb_dims for nb_dims in nb_dims_in if nb_dims > 0)
    if ind_ellipsis is not None:
        nb_dims_in[ind_ellipsis] = max(0, nb_dim - nb_known_dims)

    # transform each index into a slice
    # note that we don't need to know "stop" to update the affine matrix
    nb_ind = 0
    indices0 = indices
    indices = []
    for d, ind in enumerate(indices0):
        if isinstance(ind, slice):
            start = ind.start
            step = ind.step
            step = 1 if step is None else step
            start = 0 if (start is None and step > 0) else \
                    shape[nb_ind] - 1 if (start is None and step < 0) else \
                    shape[nb_ind] + start if start < 0 else \
                    start
            indices.append(slice(start, None, step))
            nb_ind += 1
        elif ind is Ellipsis:
            for dd in range(nb_ind, nb_ind + nb_dims_in[d]):
                start = 0
                step = 1
                indices.append(slice(start, None, step))
                nb_ind += 1
        elif is_int(ind):
            indices.append(to_int(ind))
        elif ind is None:
            assert (ind is None), "Strange index of type {}".format(type(ind))
            indices.append(None)

    # Extract shift and scale in each dimension
    shifts = []
    scales = []
    slicer = []
    shape_out = []
    for d, ind in enumerate(indices):
        # translation + scale
        if isinstance(ind, slice):
            shifts.append(ind.start)
            scales.append(ind.step)
            shape_out.append(shape[d] // abs(ind.step))
            slicer.append(slice(None))
        elif isinstance(ind, int):
            scales.append(0)
            shifts.append(ind)
            slicer.append(0)
        else:
            slicer.append(None)
            assert (ind is None), "Strange index of type {}".format(type(ind))

    # build voxel-to-voxel transformation matrix
    lin = torch.diag(torch.as_tensor(scales, **backend))
    if any(not isinstance(s, slice) for s in slicer):
        # drop/add columns
        lin = torch.unbind(lin, dim=-1)
        zero = torch.zeros(len(shifts), **backend)
        new_lin = []
        for s in slicer:
            if isinstance(s, slice):
                col, *lin = lin
                new_lin.append(col)
            elif isinstance(s, int):
                col, *lin = lin
            elif s is None:
                new_lin.append(zero)
        lin = torch.stack(new_lin, dim=-1) if new_lin else []
    trl = torch.as_tensor(shifts, **backend)[..., None]
    trf34 = torch.cat((lin, trl), dim=1) if len(lin) else trl
    trf = torch.eye(4, **backend)
    trf[:-1, :] = trf34

    # compose
    affine = affine.matmul(trf)
    return affine, tuple(shape_out)


####

def make_kernels1d():
    F0 = 151 / 315
    F1 = 397 / 1680
    F2 = 1 / 42
    F3 = 1 / 5040
    G0 = 2 / 3
    G1 = -1 / 8
    G2 = -1 / 5
    G3 = -1 / 120
    H0 = 8 / 3
    H1 = -3 / 2
    H2 = 0
    H3 = 1 / 6
    FG0 = 0
    FG1 = -49/144
    FG2 = -7/90
    FG3 = -1/720
    F = [F3, F2, F1, F0, F1, F2, F3]
    G = [G3, G2, G1, G0, G1, G2, G3]
    H = [H3, H2, H1, H0, H1, H2, H3]
    FG = [-FG3, -FG2, -FG1, FG0, FG1, FG2, FG3]
    F = torch.as_tensor(F, dtype=torch.double)
    G = torch.as_tensor(G, dtype=torch.double)
    H = torch.as_tensor(H, dtype=torch.double)
    FG = torch.as_tensor(FG, dtype=torch.double)
    return F, G, H, FG


def make_absolute3_kernel():
    F, *_ = make_kernels1d()
    K = F[None, None, :] * F[None, :, None] * F[:, None, None]
    return K


def make_membrane3_kernel():
    F, G, *_ = make_kernels1d()
    K = (F[None, None, :] * F[None, :, None] * G[:, None, None] +
         F[None, None, :] * G[None, :, None] * F[:, None, None] +
         G[None, None, :] * F[None, :, None] * F[:, None, None])
    return K


def make_bending3_kernel():
    F, G, H, *_ = make_kernels1d()
    K = (F[None, None, :] * F[None, :, None] * H[:, None, None] +
         F[None, None, :] * H[None, :, None] * F[:, None, None] +
         H[None, None, :] * F[None, :, None] * F[:, None, None] +
         F[None, None, :] * G[None, :, None] * G[:, None, None] * 2 +
         G[None, None, :] * F[None, :, None] * G[:, None, None] * 2 +
         G[None, None, :] * G[None, :, None] * F[:, None, None] * 2)
    return K


def make_linearelastic3_kernel():
    FF, GG, HH, FG = make_kernels1d()
    # diagonal of lam (divergence)
    Kxx = GG[None, None, :] * FF[None, :, None] * FF[:, None, None]
    Kyy = FF[None, None, :] * GG[None, :, None] * FF[:, None, None]
    Kzz = FF[None, None, :] * FF[None, :, None] * GG[:, None, None]
    # off diagonal (common to lam and mu)
    Kxy = - FG[None, None, :] * FG[None, :, None] * FF[:, None, None]
    Kxz = - FG[None, None, :] * FF[None, :, None] * FG[:, None, None]
    Kyz = - FF[None, None, :] * FG[None, :, None] * FG[:, None, None]
    # diagonal of mu == membrane of each component
    return Kxx, Kyy, Kzz, Kxy, Kxz, Kyz


kernels1d = make_kernels1d()
absolute3_kernel = make_absolute3_kernel()
membrane3_kernel = make_membrane3_kernel()
bending3_kernel = make_bending3_kernel()
linearelastic3_kernel = make_linearelastic3_kernel()


# The convolution kernels are separable so could be applied as a series
# of 1D convolutions. However, I've done a quick benchmark and it does
# not look beneficial.
#
# Benchmark on a [192, 192, 192, 3] field
# GPU
#   separable = True:  230 ms
#   separable = False:  70 ms
# CPU (2 x 20 cores)
#   separable = True:  600 ms
#   separable = False: 400 ms
# so better to do one 7x7x7 convolution


def absolute3(x, bound='circular'):
    """Absolute energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * membrane3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, absolute3_kernel, bound)


def membrane3(x, bound='circular'):
    """Membrane energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * membrane3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, membrane3_kernel, bound)


def bending3(x, bound='circular'):
    """Bending energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * bending3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, bending3_kernel, bound)


def linearelastic3(x, mu=0.05, lam=0.2, bound='circular'):
    """Linear-elastic energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * linearelastic3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    mu : float
        Second lame constant (penalty on shears)
    lam : float
        First lame constant (penalty on divergence)
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    Kxx, Kyy, Kzz, Kxy, Kxz, Kyz = linearelastic3_kernel
    M = membrane3_kernel.to(x)
    K = M.new_empty([3, 3, 7, 7, 7])
    K[0, 0] = lam * Kxx + mu * M
    K[1, 1] = lam * Kyy + mu * M
    K[2, 2] = lam * Kzz + mu * M
    K[0, 1] = K[1, 0] = (lam + mu) * Kxy
    K[0, 2] = K[2, 0] = (lam + mu) * Kxz
    K[1, 2] = K[2, 1] = (lam + mu) * Kyz
    y = _conv3(x, K, bound)
    return y


def _conv3(x, kernel, bound='circular'):
    kernel = kernel.to(x)
    if kernel.ndim == 5:
        # linear elastic -> (3, 3, 7, 7, 7) kernel
        x = x.movedim(-1, 0)[None]
        y = functional.pad(x, [3]*6, mode=bound)
        y = functional.conv3d(y, kernel)
        y= y[0].movedim(0, -1)
    else:
        # absolute/membrane/bending -> (7, 7, 7) kernel
        x = x[None].movedim(-1, 0)
        y = functional.pad(x, [3]*6, mode=bound)
        y = functional.conv3d(y, kernel[None, None])
        y = y.movedim(0, -1)[0]
    return y


class LabelDataset(Dataset):

    def __init__(self, fnames):
        self.fnames = fnames

    def __len__(self):
        return len(self.fnames)

    def __getitem__(self, item):
        print(item, self.fnames[item])
        prior = sp.load_npz(self.fnames[item])
        prior_indices = torch.as_tensor(prior.row)
        prior_values = torch.as_tensor(prior.data)
        return prior_indices, prior_values


def get_priors_or_posteriors(atlas_names, atlas_size, grids, atlas_resolution,
                             resampled_resolution, tissue_index,
                             n_tissues, n_labels, aff_A, aff_A_r,
                             gaussian_lhoods=None, number_of_gaussians=None, normalizer=None,
                             LUT_file=None, aff_r=None, prefetch=4, workers=2,
                             dtype=torch.float32, numpy_dtype=np.float32, device='cpu'):

    """Load the priors in efficiently by centering the bounding box on the structure.

    The function can be used to load and lump the priors to "super classes" and to
    compute the posteriors if the gaussian likelihoods are defined.

    Parameters
    ----------
    atlas_names : list[str]
        list of the paths to each individual atlas class (.npz files)
    aff_A : (4,4) array
        Voxel-to-world transformation of the atlas.
    aff_A_r : (4,4) array
        Voxel-to-world transformation of the resampled atlas.
    aff_r : (4,4) array
        Voxel-to-world transformation of the resampled image.
    atlas_size : (3) array
        Size of the atlas
    grids : (nx, ny, nz, 3) array
        Sampling grid of the atlas in voxel space
    atlas_resolution : float
        Resolution of the atlas (assumed isotropic)
    resampled_resolution : float
        Resolution of the resampled atlas (assumed isotropic)
    tissue_index : list[int]
        Mapping of the atlas structures to tissues
    hemi: str
        String denoting the hemisphere "l" or "r" (default "l")
    gaussian_lhoods : (nx, ny, nz, num_gaussian) array
        Gaussian likelihoods of all mixture components. num_gaussians
        is sum(number_of_gaussians)
    number_of_gaussians : list[int]
        Number of mixture components for each tissue class
    normalizer : (nx, ny, nz) array
        Normalizer for the posterior
    n_labels: int
        Number of labels in the atlas
    prefetch : int
        How many priors to preload using the LabelDataset class
    workers : int
         How many workers to use to process the priors


    """
    # TODO: choose good number of workers/prefetch factor
    if gaussian_lhoods is not None:
        seg = torch.zeros(normalizer.shape, dtype=torch.int, device=device)
        seg_rgb = torch.zeros([*normalizer.shape, 3], dtype=dtype, device=device)
        max_p = torch.zeros(normalizer.shape, dtype=dtype, device=device)
        vols = torch.zeros(n_labels, device=device, dtype=dtype)
        names, colors = my.read_LUT(LUT_file)
        voxel_vol = np.abs(np.linalg.det(aff_r))
    else:
        A = np.zeros([*atlas_size, n_tissues], dtype=numpy_dtype)

    resolutions_match = resampled_resolution == atlas_resolution
    prefetch_factor = max(prefetch//workers, 1)
    label_loader = DataLoader(LabelDataset(atlas_names), num_workers=workers, prefetch_factor=prefetch_factor)
    for n, (prior_indices, prior_values) in enumerate(label_loader):
        print('Reading in label ' + str(n + 1) + ' of ' + str(n_labels), end='\r', flush=True)

        if prior_indices.numel() == 0:
            continue
        prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
        prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()

        aff_pr = np.copy(aff_A)
        if n == 0:
            # background
            prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                            [torch.Size(atlas_size).numel()]).to_dense()
            del prior_indices, prior_values
            prior = prior.reshape(torch.Size(atlas_size))
            lr_crop = (slice(None),) * 3
        else:
            # find bounding box of label in atlas space
            prior_indices = ind2sub(prior_indices, atlas_size)
            min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item() + 1
            min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item() + 1
            min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item() + 1
            crop_atlas_size = [max_x - min_x, max_y - min_y, max_z - min_z]
            prior_indices[0] -= min_x
            prior_indices[1] -= min_y
            prior_indices[2] -= min_z
            prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
            del prior_indices, prior_values
            aff_pr[:3, -1] += aff_pr[:3, :3] @ np.asarray([min_x, min_y, min_z])
            # find bounding box of label in MRI space
            hr2lr = np.linalg.inv(aff_A_r) @ aff_A
            min_x, min_y, min_z = (hr2lr[:3, :3] @ np.asarray([min_x-1, min_y-1, min_z-1] + hr2lr[:3, -1])).tolist()
            max_x, max_y, max_z = (hr2lr[:3, :3] @ np.asarray([max_x, max_y, max_z] + hr2lr[:3, -1])).tolist()
            mask =  (grids[0, ..., 0] >= min_x)
            mask &= (grids[0, ..., 0] <= max_x)
            mask &= (grids[0, ..., 1] >= min_y)
            mask &= (grids[0, ..., 1] <= max_y)
            mask &= (grids[0, ..., 2] >= min_z)
            mask &= (grids[0, ..., 2] <= max_z)
            if ~mask.any():
                continue
            nx, ny, nz = mask.shape
            tmp = mask.reshape([nx, -1]).any(-1).nonzero()
            lr_min_x, lr_max_x = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.movedim(0, -1).reshape([ny, -1]).any(-1).nonzero()
            lr_min_y, lr_max_y = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.reshape([-1, nz]).any(0).nonzero()
            lr_min_z, lr_max_z = tmp.min().item(), tmp.max().item() + 1
            del tmp, mask
            lr_crop = (slice(lr_min_x, lr_max_x), slice(lr_min_y, lr_max_y), slice(lr_min_z, lr_max_z))

        if not resolutions_match:
            prior, aff_pr = torch_resize(prior, aff_pr, resampled_resolution, device, dtype=dtype)

        # shift/scale sampling grid appropriately
        if not (resolutions_match and n==0):
            aff_shift = torch.as_tensor(np.linalg.inv(aff_pr) @ aff_A_r, device=device, dtype=dtype)
            grids_shifted = grids[(slice(None), *lr_crop, slice(None))]
            grids_shifted = aff_shift[:3, :3].matmul(grids_shifted.unsqueeze(-1)).squeeze(-1)
            grids_shifted = grids_shifted.add_(aff_shift[:3, -1])
            breakpoint()
            prior = interpol.grid_pull(prior[None, None, :, :, :], grids_shifted, interpolation=1)
            del grids_shifted

        if gaussian_lhoods is not None:
            num_gaussians = number_of_gaussians[tissue_index[n]]

            num_components = number_of_gmm_components[c+1]
            gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[tissue_index[n]]) + \
                                            np.array(range(num_components)), device=device, dtype=dtype).int()
            lhood = torch.sum(gaussian_lhoods[:, :, :, gaussian_numbers], 3)
            post = torch.squeeze(prior)

            post *= lhood[lr_crop]
            post /= normalizer[lr_crop]
            del prior

            vols[n] = torch.sum(post) * voxel_vol
            mask = (post > max_p[lr_crop])
            max_p[lr_crop][mask] = post[mask]
            lab = int(label_list[n])
            seg[lr_crop].masked_fill_(mask, lab)
            del mask
            for c in range(3):
                seg_rgb[(*lr_crop, c)].add_(post, alpha=colors[lab][c])
        else:
            A[(*lr_crop, tissue_index[n])] = A[(*lr_crop, tissue_index[n])] + prior.cpu().numpy()

    if gaussian_lhoods is not None:
        return seg, seg_rgb, vols
    else:
        return A
