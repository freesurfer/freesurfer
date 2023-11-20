import nibabel as nib
import torch
import numpy as np
import os
from scipy.ndimage import gaussian_filter

#######

def viewVolume(x, aff=None):

    if aff is None:
        aff = np.eye(4)
    else:
        if type(aff) == torch.Tensor:
            aff = aff.cpu().detach().numpy()

    if type(x) is not list:
        x = [x]

    cmd = 'source /usr/local/freesurfer/nmr-dev-env-bash && freeview '

    for n in range(len(x)):
        vol = x[n]
        if type(vol) == torch.Tensor:
            vol = vol.cpu().detach().numpy()
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

###############################

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

###############################

def myzoom_torch(X, factor, device, aff=None):

    if len(X.shape)==3:
        X = X[..., None]

    delta = (1.0 - factor) / (2.0 * factor)
    newsize = np.round(X.shape[:-1] * factor).astype(int)

    vx = torch.arange(delta[0], delta[0] + newsize[0] / factor[0], 1 / factor[0], dtype=torch.float, device=device)[:newsize[0]]
    vy = torch.arange(delta[1], delta[1] + newsize[1] / factor[1], 1 / factor[1], dtype=torch.float, device=device)[:newsize[1]]
    vz = torch.arange(delta[2], delta[2] + newsize[2] / factor[2], 1 / factor[2], dtype=torch.float, device=device)[:newsize[2]]

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

    for channel in range(X.shape[3]):
        Xc = X[:,:,:,channel]

        tmp1 = torch.zeros([newsize[0], Xc.shape[1], Xc.shape[2]], dtype=torch.float, device=device)
        for i in range(newsize[0]):
            tmp1[i, :, :] = wfx[i] * Xc[fx[i], :, :] +  wcx[i] * Xc[cx[i], :, :]
        tmp2 = torch.zeros([newsize[0], newsize[1], Xc.shape[2]], dtype=torch.float, device=device)
        for j in range(newsize[1]):
            tmp2[:, j, :] = wfy[j] * tmp1[:, fy[j], :] +  wcy[j] * tmp1[:, cy[j], :]
        for k in range(newsize[2]):
            Y[:, :, k, channel] = wfz[k] * tmp2[:, :, fz[k]] +  wcz[k] * tmp2[:, :, cz[k]]

    if Y.shape[3] == 1:
        Y = Y[:,:,:, 0]

    if aff is not None:
        aff_new = aff.copy()
        for c in range(3):
            aff_new[:-1, c] = aff_new[:-1, c] / factor
        aff_new[:-1, -1] = aff_new[:-1, -1] - aff[:-1, :-1] @ (0.5 - 0.5 / (factor * np.ones(3)))
        return Y, aff_new
    else:
        return Y

############################


def myzoom_torch_anisotropic(X, aff, newsize, device):

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

############################

def torch_resize(I, aff, resolution, device, power_factor_at_half_width=5, dtype=torch.float32, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_resize(I, aff, resolution, device, power_factor_at_half_width, dtype, slow)

    slow = slow or (device == 'cpu')
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
            It = It.permute([0, 1, 3, 4, 2])
            if sigmas[d]>0:
                sl = np.ceil(sigmas[d] * 2.5).astype(int)
                v = np.arange(-sl, sl + 1)
                gauss = np.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / np.sum(gauss)
                kernel = torch.tensor(kernel,  device=device, dtype=dtype)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None, None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])


        It = torch.squeeze(It)
        It, aff2 = myzoom_torch_anisotropic(It, aff, newsize, device)
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



#######


def align_volume_to_ref(volume, aff, aff_ref=None, return_aff=False, n_dims=3):
    """This function aligns a volume to a reference orientation (axis and direction) specified by an affine matrix.
    :param volume: a numpy array
    :param aff: affine matrix of the floating volume
    :param aff_ref: (optional) affine matrix of the target orientation. Default is identity matrix.
    :param return_aff: (optional) whether to return the affine matrix of the aligned volume
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: aligned volume, with corresponding affine matrix if return_aff is True.
    """

    # work on copy
    aff_flo = aff.copy()

    # default value for aff_ref
    if aff_ref is None:
        aff_ref = np.eye(4)

    # extract ras axes
    ras_axes_ref = get_ras_axes(aff_ref, n_dims=n_dims)
    ras_axes_flo = get_ras_axes(aff_flo, n_dims=n_dims)

    # align axes
    aff_flo[:, ras_axes_ref] = aff_flo[:, ras_axes_flo]
    for i in range(n_dims):
        if ras_axes_flo[i] != ras_axes_ref[i]:
            volume = torch.swapaxes(volume, ras_axes_flo[i], ras_axes_ref[i])
            swapped_axis_idx = np.where(ras_axes_flo == ras_axes_ref[i])
            ras_axes_flo[swapped_axis_idx], ras_axes_flo[i] = ras_axes_flo[i], ras_axes_flo[swapped_axis_idx]

    # align directions
    dot_products = np.sum(aff_flo[:3, :3] * aff_ref[:3, :3], axis=0)
    for i in range(n_dims):
        if dot_products[i] < 0:
            # volume = torch.from_numpy(volume)
            volume = torch.flip(volume, [i])
            aff_flo[:, i] = - aff_flo[:, i]
            aff_flo[:3, 3] = aff_flo[:3, 3] - aff_flo[:3, i] * (volume.shape[i] - 1)

    if return_aff:
        return volume, aff_flo
    else:
        return volume

##############

def get_ras_axes(aff, n_dims=3):
    """This function finds the RAS axes corresponding to each dimension of a volume, based on its affine matrix.
    :param aff: affine matrix Can be a 2d numpy array of size n_dims*n_dims, n_dims+1*n_dims+1, or n_dims*n_dims+1.
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: two numpy 1d arrays of lengtn n_dims, one with the axes corresponding to RAS orientations,
    and one with their corresponding direction.
    """
    aff_inverted = np.linalg.inv(aff)
    img_ras_axes = np.argmax(np.absolute(aff_inverted[0:n_dims, 0:n_dims]), axis=0)
    return img_ras_axes
