import os
import numpy as np
import scipy.io as sio
import prep

def load_surf_by_name(surf_name, curv='sulc')

    (path_name, name) = os.path.split(surf_name)
    hemi = name[0:3]
    curv_name = path_name + '/' + hemi + curv
    X = spherePad(surf_name, curv_name)
    X = np.reshape(X, (1,) + X.shape + (1,))

    return X


def sphere_gen(surf_names):
    #idx = 0
    while(True):
        idx = np.random.randint(len(surf_names))
        X = load_surf_by_name(surf_names[idx])
        return_vals = [X]

        yield tuple(return_vals)


def load_example_by_name(vol_name, seg_name):

    X = np.load(vol_name)['vol_data']
    X = np.reshape(X, (1,) + X.shape + (1,))

    return_vals = [X]

    X_seg = np.load(seg_name)['vol_data']
    X_seg = np.reshape(X_seg, (1,) + X_seg.shape + (1,))
    return_vals.append(X_seg)

    return tuple(return_vals)


def load_img_by_name(vol_name, seg_name=None):

    X = sio.loadmat(vol_name)['tmp']
    X = np.reshape(X, (1,) + X.shape + (1,))

    #return_vals = [X]

    #X_seg = np.load(seg_name)['vol_data']
    #X_seg = np.reshape(X_seg, (1,) + X_seg.shape + (1,))
    #return_vals.append(X_seg)

    return X


def example_gen(vol_names, return_segs=False, seg_dir=None):
    #idx = 0
    while(True):
        idx = np.random.randint(len(vol_names))
        X = np.load(vol_names[idx])['vol_data']
        X = np.reshape(X, (1,) + X.shape + (1,))

        return_vals = [X]

        if(return_segs):
            name = os.path.basename(vol_names[idx])
            X_seg = np.load(seg_dir + name[0:-8]+'aseg.npz')['vol_data']
            X_seg = np.reshape(X_seg, (1,) + X_seg.shape + (1,))
            return_vals.append(X_seg)

        # print vol_names[idx] + "," + seg_dir + name[0:-8]+'aseg.npz'

        yield tuple(return_vals)

def sphericaldata_gen(vol_names, return_segs=False, seg_dir=None):
    #idx = 0
    while(True):
        idx = np.random.randint(len(vol_names))
        X = sio.loadmat(vol_names[idx])['tmp']
        X = np.reshape(X, (1,) + X.shape + (1,))

        return_vals = [X]

        if(return_segs):
            name = os.path.basename(vol_names[idx])
            X_seg = np.load(seg_dir + name[0:-8]+'aseg.npz')['vol_data']
            X_seg = np.reshape(X_seg, (1,) + X_seg.shape + (1,))
            return_vals.append(X_seg)

        # print vol_names[idx] + "," + seg_dir + name[0:-8]+'aseg.npz'

        yield tuple(return_vals)
