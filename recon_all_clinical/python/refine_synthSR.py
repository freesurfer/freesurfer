# python imports
import os
import sys
import numpy as np
from argparse import ArgumentParser
import nibabel as nib
from scipy.ndimage import binary_dilation
from scipy.ndimage.morphology import distance_transform_edt

# add main folder to python path and import SynthSR packages
code_home = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
print(code_home)
sys.path.append(code_home)




# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================

def main():

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("input", type=str,
                        help="SynthSR input")
    parser.add_argument("ribbon", type=str,
                        help="Ribbon file")
    parser.add_argument("fake", type=str,
                        help="fake image file")
    parser.add_argument("output", type=str,
                        help="output")

    args = vars(parser.parse_args())

    print('Reading inputs')
    im, aff, hdr = load_volume(args['input'], im_only=False, dtype='float')
    imR, affR, hdrR = load_volume(args['ribbon'], im_only=False, dtype='float')
    imF, affF, hdrF = load_volume(args['fake'], im_only=False, dtype='float')


    if np.max(np.abs(aff-affR))>1e-3:
        raise Exception('headers differ')
    if np.max(np.abs(aff-affF))>1e-3:
        raise Exception('headers differ')

    print('Processing')
    M = ( (imR==3) | (imR==42) )
    D = distance_transform_edt(M==0)
    D = D - 1.0
    D[D<0] = 0
    Prib = np.exp(-1.0*D)
    Pfake = (1.0 - Prib)
    I = (imF * Prib) + (im * Pfake)

    print('Writing to disk')
    save_volume(I, aff, None, args['output'])

    print('Done!')



# ================================================================================================
#                                         Auxiliary functions
# ================================================================================================

def save_volume(volume, aff, header, path, res=None, dtype=None, n_dims=3):

    mkdir(os.path.dirname(path))
    if '.npz' in path:
        np.savez_compressed(path, vol_data=volume)
    else:
        if header is None:
            header = nib.Nifti1Header()
        if isinstance(aff, str):
            if aff == 'FS':
                aff = np.array([[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]])
        elif aff is None:
            aff = np.eye(4)
        nifty = nib.Nifti1Image(volume, aff, header)
        if dtype is not None:
            if 'int' in dtype:
                volume = np.round(volume)
            volume = volume.astype(dtype=dtype)
            nifty.set_data_dtype(dtype)
        if res is not None:
            if n_dims is None:
                n_dims, _ = get_dims(volume.shape)
            res = reformat_to_list(res, length=n_dims, dtype=None)
            nifty.header.set_zooms(res)
        nib.save(nifty, path)


def mkdir(path_dir):
    """Recursively creates the current dir as well as its parent folders if they do not already exist."""
    if path_dir[-1] == '/':
        path_dir = path_dir[:-1]
    if not os.path.isdir(path_dir):
        list_dir_to_create = [path_dir]
        while not os.path.isdir(os.path.dirname(list_dir_to_create[-1])):
            list_dir_to_create.append(os.path.dirname(list_dir_to_create[-1]))
        for dir_to_create in reversed(list_dir_to_create):
            os.mkdir(dir_to_create)

def build_binary_structure(connectivity, n_dims, shape=None):
    """Return a dilation/erosion element with provided connectivity"""
    if shape is None:
        shape = [connectivity * 2 + 1] * n_dims
    else:
        shape = reformat_to_list(shape, length=n_dims)
    dist = np.ones(shape)
    center = tuple([tuple([int(s / 2)]) for s in shape])
    dist[center] = 0
    dist = distance_transform_edt(dist)
    struct = (dist <= connectivity) * 1
    return struct

def load_volume(path_volume, im_only=True, squeeze=True, dtype=None, aff_ref=None):

    assert path_volume.endswith(('.nii', '.nii.gz', '.mgz', '.npz')), 'Unknown data file: %s' % path_volume

    if path_volume.endswith(('.nii', '.nii.gz', '.mgz')):
        x = nib.load(path_volume)
        if squeeze:
            volume = np.squeeze(x.get_fdata())
        else:
            volume = x.get_fdata()
        aff = x.affine
        header = x.header
    else:  # npz
        volume = np.load(path_volume)['vol_data']
        if squeeze:
            volume = np.squeeze(volume)
        aff = np.eye(4)
        header = nib.Nifti1Header()
    if dtype is not None:
        if 'int' in dtype:
            volume = np.round(volume)
        volume = volume.astype(dtype=dtype)

    # align image to reference affine matrix
    if aff_ref is not None:
        n_dims, _ = get_dims(list(volume.shape), max_channels=10)
        volume, aff = align_volume_to_ref(volume, aff, aff_ref=aff_ref, return_aff=True, n_dims=n_dims)

    if im_only:
        return volume
    else:
        return volume, aff, header

# execute script
if __name__ == '__main__':
    main()
