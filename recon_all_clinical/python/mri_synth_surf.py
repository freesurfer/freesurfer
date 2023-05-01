
# python imports
import os
import sys
import numpy as np
from argparse import ArgumentParser
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend
backend.set_image_data_format('channels_last')
import nibabel as nib
from scipy.ndimage.morphology import distance_transform_edt
from scipy.ndimage import binary_dilation
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

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
    parser.add_argument("--input_image", type=str,
                        help="images to process. Can be the path to a single image or to a folder")
    parser.add_argument("--input_synthseg", type=str,
                        help="path to synthseg segmentation")
    parser.add_argument("--subject_mri_dir", type=str,
                        help="path of mri directory of this subject, e.g., $FREESURFER_HOME/subjects/bert/mri/")
    parser.add_argument("--model_file", type=str,
                        help="path to model file")
    parser.add_argument("--cpu", action="store_true", help="enforce running with CPU rather than GPU.")
    parser.add_argument("--threads", type=int, default=1, dest="threads",
                        help="number of threads to be used by tensorflow when running on CPU.")
    parser.add_argument("--pad", type=int, default=0, dest="pad",
                        help="pad with these many voxels in every direction")

    args = vars(parser.parse_args())

    # enforce CPU processing if necessary
    if args['cpu']:
        print('using CPU, hiding all CUDA_VISIBLE_DEVICES')
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

    # limit the number of threads to be used if running on CPU
    tf.config.threading.set_intra_op_parallelism_threads(args['threads'])

    # mri dir
    mri_dir = os.path.abspath(args['subject_mri_dir'])

    print('Reading image')
    im, aff, hdr = load_volume(args['input_image'], im_only=False, dtype='float')

    print('Reading SynthSeg segmentation / parcellation')
    imS, affS, hdrS = load_volume(args['input_synthseg'], im_only=False, dtype='float')

    print('Resampling and padding image')
    while len(im.shape)>3: # in case it's rgb
        im = np.mean(im, axis=-1)
    im, aff = resample_volume(im, aff, [1.0, 1.0, 1.0])
    im, aff2 = align_volume_to_ref(im, aff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
    im = im - np.min(im)
    im = im / np.max(im)
    I = im[np.newaxis, ..., np.newaxis]
    W = (np.ceil(np.array(I.shape[1:-1]) / 32.0) * 32).astype('int')
    idx = np.floor((W - I.shape[1:-1]) / 2).astype('int')
    S = np.zeros([1, *W, 1])
    S[0, idx[0]:idx[0] + I.shape[1], idx[1]:idx[1] + I.shape[2], idx[2]:idx[2] + I.shape[3], :] = I

    print('Building U-net')
    unet_model = unet(nb_features=24,
                     input_shape=[None, None, None, 1],
                     nb_levels=5,
                     conv_size=3,
                     nb_labels=9,
                     feat_mult=2,
                     nb_conv_per_level=2,
                     conv_dropout=0,
                     final_pred_activation='linear',
                     batch_norm=-1,
                     activation='elu',
                     input_model=None)

    unet_model.load_weights(args['model_file'], by_name=True)

    print('Pushing image through U-net')
    output = unet_model.predict(S)
    pred = np.squeeze(output)
    pred = pred[idx[0]:idx[0] + I.shape[1], idx[1]:idx[1] + I.shape[2], idx[2]:idx[2] + I.shape[3], :]

    print('Saving prediction to disk')
    prediction_file = mri_dir + '/synthsurf.mgz'
    save_volume(pred, aff2, None, prediction_file)

    print('Resampling SynthSeg to space of prediction')
    II, JJ, KK = np.meshgrid(np.arange(pred.shape[0]), np.arange(pred.shape[1]), np.arange(pred.shape[2]), indexing='ij')
    M = np.matmul(np.linalg.inv(affS), aff2)
    II2 = M[0, 0] * II + M[0, 1] * JJ + M[0, 2] * KK + M[0, 3]
    JJ2 = M[1, 0] * II + M[1, 1] * JJ + M[1, 2] * KK + M[1, 3]
    KK2 = M[2, 0] * II + M[2, 1] * JJ + M[2, 2] * KK + M[2, 3]
    IIr = np.round(II2).astype(int)
    JJr = np.round(JJ2).astype(int)
    KKr = np.round(KK2).astype(int)
    IIr[IIr < 0] = 0
    JJr[JJr < 0] = 0
    KKr[KKr < 0] = 0
    IIr[IIr > (imS.shape[0] - 1)] = (imS.shape[0] - 1)
    JJr[JJr > (imS.shape[1] - 1)] = (imS.shape[1] - 1)
    KKr[KKr > (imS.shape[2] - 1)] = (imS.shape[2] - 1)
    S = imS[IIr, JJr, KKr]

    print('Computing Talairach transform')
    # Note: labels and centers of gravity precomputed from:
    # /autofs/space/panamint_005/users/iglesias/data/MNItemplates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c.nii.gz
    # /autofs/space/panamint_005/users/iglesias/data/MNItemplates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c.synthseg2.nii.gz
    labels=np.array([ 2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26, 28, 41, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 1001, 1002, 1003, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035]).astype(int)
    atlCOG = np.array([[-27., -13., -31., -17., -24., -11., -14., -27., -20.,   0.,   0.,
              0., -26., -23.,   0.,  -8.,  -9.,  27.,  14.,  33.,  17.,  25.,
             11.,  14.,  27.,  20.,  26.,  23.,   9.,   9., -54.,  -4., -37.,
             -6., -23., -35., -44., -53.,  -6., -30., -23., -13.,  -5., -60.,
            -23.,  -6., -49., -44., -49., -11., -48.,  -4., -44.,  -8.,  -4.,
            -34., -10., -25., -54., -57.,  -9., -28., -45., -37.,  55.,   4.,
             37.,   6.,  22.,  35.,  45.,  52.,   6.,  31.,  23.,  13.,   5.,
             60.,  23.,   6.,  49.,  43.,  50.,  10.,  48.,   4.,  43.,   8.,
              4.,  33.,  10.,  25.,  54.,  57.,   8.,  30.,  46.,  37.],
           [-18., -18., -13., -54., -63., -18.,  12.,   3.,  -2.,  -7., -46.,
            -31., -20.,  -4., -21.,  10., -16., -18., -21., -15., -54., -63.,
            -18.,  12.,   3.,  -3., -20.,  -4.,  10., -16., -43.,  21.,  11.,
            -79.,  -4., -41., -67., -35., -46., -90.,  29., -67.,  42., -23.,
            -32., -28.,  16.,  43.,  32., -79., -19., -21.,  -7., -59.,  38.,
             49.,  29., -63.,  -9., -35.,  67.,  13., -21.,   2., -42.,  20.,
             13., -81.,  -6., -41., -66., -34., -46., -90.,  29., -68.,  41.,
            -24., -32., -27.,  16.,  44.,  32., -82., -19., -20.,  -7., -59.,
             39.,  50.,  30., -62.,  -8., -34.,  67.,  13., -20.,   2.],
           [ 19.,  15., -15., -35., -38.,   6.,  10.,  -1.,  -2.,  -4., -34.,
            -34., -16., -20.,   8.,  -9., -10.,  19.,  14., -14., -35., -38.,
              6.,  10.,  -1.,  -3., -16., -20.,  -9., -10.,   8.,  28.,  49.,
             20., -34., -21.,  31., -24.,  22.,   0., -19.,  -5., -17., -13.,
            -17.,  58.,  14., -14.,   2.,   7.,  47.,  40.,  47.,  38.,   1.,
             20.,  47.,  53.,  -4.,  33., -11., -37.,   9.,  -2.,   6.,  28.,
             48.,  21., -33., -21.,  31., -25.,  23.,  -1., -20.,  -5., -17.,
            -13., -17.,  59.,  13., -14.,   2.,   7.,  46.,  40.,  47.,  39.,
            1.,  19.,  47.,  54.,  -5.,  34., -11., -38.,   8.,  -3.]]).astype(float)
    nlab = len(labels)

    refCOG = np.zeros([4, nlab])
    ok = np.ones(nlab)
    for l in range(nlab):
        aux = np.where(imS == labels[l])
        if len(aux[0]) > 50:
            refCOG[0, l] = np.median(aux[0])
            refCOG[1, l] = np.median(aux[1])
            refCOG[2, l] = np.median(aux[2])
            refCOG[3, l] = 1
        else:
            ok[l] = 0
    refCOG = np.matmul(affS, refCOG)[:-1, :]  # in RAS
    TAL = getM(refCOG[:,ok>0], atlCOG[:,ok>0])
    with open(mri_dir + '/transforms/talairach.xfm', 'w') as f:
        f.write('MNI Transform File')
        f.write('\n')
        f.write('% avi2talxfm')
        f.write('\n')
        f.write('\n')
        f.write('Transform_Type = Linear;')
        f.write('\n')
        f.write('Linear_Transform = ')
        f.write('\n')
        f.write(str(TAL[0,0]) + ' ' + str(TAL[0,1]) + ' ' + str(TAL[0,2]) + ' ' + str(TAL[0,3]))
        f.write('\n')
        f.write(str(TAL[1,0]) + ' ' + str(TAL[1,1]) + ' ' + str(TAL[1,2]) + ' ' + str(TAL[1,3]))
        f.write('\n')
        f.write(str(TAL[2,0]) + ' ' + str(TAL[2,1]) + ' ' + str(TAL[2,2]) + ' ' + str(TAL[2,3]) + ';')
        f.write('\n')

    print('Making synthetic image and other <<fake>> volumes, in FreeSurfer orientation!')
    I = im.copy()
    # Distance parameter
    a = 2
    # % Create two fake images, one per hemi
    W = pred[:,:,:,0]
    P = pred[:,:,:,1]
    # Fleft = 70*(1-(np.tanh(a*W)+1)/2) + 40*(1-(np.tanh(a*P)+1)/2) TODO: undo this hack? It's for cosmetic purposes...
    Fleft = 70 * (1 - (np.tanh(a * (W+0.3)) + 1) / 2) + 40 * (1 - (np.tanh(a * P) + 1) / 2)
    W = pred[:,:,:,2]
    P = pred[:,:,:,3]
    # Fright = 70*(1-(np.tanh(a*W)+1)/2) + 40*(1-(np.tanh(a*P)+1)/2) TODO: undo this hack? It's for cosmetic purposes...
    Fright = 70 * (1 - (np.tanh(a * (W+0.3)) + 1) / 2) + 40 * (1 - (np.tanh(a * P) + 1) / 2)
    # Split whole volume into left and right using distance maps, and use the
    # division to compose the two fake images into one
    Scopy = S.copy()
    S[(S==14) | (S==15) | (S==16) | (S==24)] = 0 # kill non-lateral labels
    L = ((S>0) &  (S<30)) | ((S>1000) & (S<2000))
    R = (S>2000) | ((S<1000) & (S>30))
    Dleft = distance_transform_edt(L==0)
    Dright = distance_transform_edt(R==0)
    L = Dleft<Dright
    R = (L==0)
    F = np.zeros(R.shape)
    F[L] = Fleft[L]
    F[R] = Fright[R]
    # kill cerebellum, binarize segmentation, dilate by 3, and kill stuff outside
    S[(S==7) | (S==8) | (S==8) | (S==46) | (S==47)] = 0
    dilate = 3
    dim = 3
    M = binary_dilation((S>0), build_binary_structure(dilate, dim))
    F[M==0] = 0

    # We also padd a bit when saving...
    aff_FS = np.array([[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]).astype(float)
    imSave, affSave = align_volume_to_ref(F, aff2, aff_ref=aff_FS, return_aff=True, n_dims=3)

    pad = args['pad']
    imSavePad = np.zeros((np.array(imSave.shape) + 2*pad).astype(int))
    imSavePad[pad:-pad, pad:-pad, pad:-pad] = imSave
    affSavePad = affSave.copy()
    affSavePad[:3, -1] = affSavePad[:3, -1] - (affSavePad[:3, :-1]  @ np.array([pad, pad, pad]))
    save_volume(imSavePad, affSavePad, None, mri_dir + '/norm.mgz')
    save_volume(imSavePad, affSavePad, None, mri_dir + '/brain.mgz')
    save_volume(imSavePad, affSavePad, None, mri_dir + '/brainmask.mgz')

    aseg = Scopy.copy()
    aseg[aseg > 2000] = 42
    aseg[aseg > 1000] = 3
    imSave, affSave = align_volume_to_ref(aseg, aff2, aff_ref=aff_FS, return_aff=True, n_dims=3)
    imSavePad[pad:-pad, pad:-pad, pad:-pad] = imSave
    save_volume(imSavePad, affSavePad, None, mri_dir + '/aseg.auto_noCCseg.mgz')

    S[S > 2000] = 42
    S[S > 1000] = 3
    WM = S.copy()
    WM[S==3] = 0
    WM[S==42] = 0
    WM[WM>0] = 110
    WM[(S==4) | (S==5) | (S==43) | (S==44)] = 250
    imSave, affSave = align_volume_to_ref(WM, aff2, aff_ref=aff_FS, return_aff=True, n_dims=3)
    imSavePad[pad:-pad, pad:-pad, pad:-pad] = imSave
    save_volume(imSavePad, affSavePad, None, mri_dir + '/wm.seg.mgz')

    FILLED = WM.copy()
    FILLED[FILLED>0] = 127
    FILLED[(FILLED>0) & (Dleft<Dright)] = 255
    imSave, affSave = align_volume_to_ref(FILLED, aff2, aff_ref=aff_FS, return_aff=True, n_dims=3)
    imSavePad[pad:-pad, pad:-pad, pad:-pad] = imSave
    save_volume(imSavePad, affSavePad, None, mri_dir + '/filled.mgz')

    print(' ')
    print('All done!')
    print(' ')



# ================================================================================================
#                                         Auxiliary functions
# ================================================================================================

# Auxiliary function to fit affine matrices
def getM(ref, mov):

   zmat = np.zeros(ref.shape[::-1])
   zcol = np.zeros([ref.shape[1], 1])
   ocol = np.ones([ref.shape[1], 1])
   zero = np.zeros(zmat.shape)

   A = np.concatenate([
       np.concatenate([np.transpose(ref), zero, zero, ocol, zcol, zcol], axis=1),
       np.concatenate([zero, np.transpose(ref), zero, zcol, ocol, zcol], axis=1),
       np.concatenate([zero, zero, np.transpose(ref), zcol, zcol, ocol], axis=1)], axis=0)

   b = np.concatenate([np.transpose(mov[0, :]), np.transpose(mov[1, :]), np.transpose(mov[2, :])], axis=0)

   x = np.matmul(np.linalg.inv(np.matmul(np.transpose(A), A)), np.matmul(np.transpose(A), b))

   M = np.stack([
       [x[0], x[1], x[2], x[9]],
       [x[3], x[4], x[5], x[10]],
       [x[6], x[7], x[8], x[11]],
       [0, 0, 0, 1]])

   return M


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

def resample_volume(volume, aff, new_vox_size, interpolation='linear'):
    """This function resizes the voxels of a volume to a new provided size, while adjusting the header to keep the RAS
    :param volume: a numpy array
    :param aff: affine matrix of the volume
    :param new_vox_size: new voxel size (3 - element numpy vector) in mm
    :return: new volume and affine matrix
    """

    pixdim = np.sqrt(np.sum(aff * aff, axis=0))[:-1]
    new_vox_size = np.array(new_vox_size)
    factor = pixdim / new_vox_size
    sigmas = 0.25 / factor
    sigmas[factor > 1] = 0  # don't blur if upsampling

    volume_filt = gaussian_filter(volume, sigmas)

    # volume2 = zoom(volume_filt, factor, order=1, mode='reflect', prefilter=False)
    x = np.arange(0, volume_filt.shape[0])
    y = np.arange(0, volume_filt.shape[1])
    z = np.arange(0, volume_filt.shape[2])

    my_interpolating_function = RegularGridInterpolator((x, y, z), volume_filt, method=interpolation)

    start = - (factor - 1) / (2 * factor)
    step = 1.0 / factor
    stop = start + step * np.ceil(volume_filt.shape * factor)

    xi = np.arange(start=start[0], stop=stop[0], step=step[0])
    yi = np.arange(start=start[1], stop=stop[1], step=step[1])
    zi = np.arange(start=start[2], stop=stop[2], step=step[2])
    xi[xi < 0] = 0
    yi[yi < 0] = 0
    zi[zi < 0] = 0
    xi[xi > (volume_filt.shape[0] - 1)] = volume_filt.shape[0] - 1
    yi[yi > (volume_filt.shape[1] - 1)] = volume_filt.shape[1] - 1
    zi[zi > (volume_filt.shape[2] - 1)] = volume_filt.shape[2] - 1

    xig, yig, zig = np.meshgrid(xi, yi, zi, indexing='ij', sparse=True)
    volume2 = my_interpolating_function((xig, yig, zig))

    aff2 = aff.copy()
    for c in range(3):
        aff2[:-1, c] = aff2[:-1, c] / factor[c]
    aff2[:-1, -1] = aff2[:-1, -1] - np.matmul(aff2[:-1, :-1], 0.5 * (factor - 1))

    return volume2, aff2

def align_volume_to_ref(volume, aff, aff_ref=None, return_aff=False, n_dims=None, return_copy=True):
    """This function aligns a volume to a reference orientation (axis and direction) specified by an affine matrix.
    :param volume: a numpy array
    :param aff: affine matrix of the floating volume
    :param aff_ref: (optional) affine matrix of the target orientation. Default is identity matrix.
    :param return_aff: (optional) whether to return the affine matrix of the aligned volume
    :param n_dims: (optional) number of dimensions (excluding channels) of the volume. If not provided, n_dims will be
    inferred from the input volume.
    :return: aligned volume, with corresponding affine matrix if return_aff is True.
    """

    # work on copy
    new_volume = volume.copy() if return_copy else volume
    aff_flo = aff.copy()

    # default value for aff_ref
    if aff_ref is None:
        aff_ref = np.eye(4)

    # extract ras axes
    if n_dims is None:
        n_dims, _ = get_dims(new_volume.shape)
    ras_axes_ref = get_ras_axes(aff_ref, n_dims=n_dims)
    ras_axes_flo = get_ras_axes(aff_flo, n_dims=n_dims)

    # align axes
    aff_flo[:, ras_axes_ref] = aff_flo[:, ras_axes_flo]
    for i in range(n_dims):
        if ras_axes_flo[i] != ras_axes_ref[i]:
            new_volume = np.swapaxes(new_volume, ras_axes_flo[i], ras_axes_ref[i])
            swapped_axis_idx = np.where(ras_axes_flo == ras_axes_ref[i])
            ras_axes_flo[swapped_axis_idx], ras_axes_flo[i] = ras_axes_flo[i], ras_axes_flo[swapped_axis_idx]

    # align directions
    dot_products = np.sum(aff_flo[:3, :3] * aff_ref[:3, :3], axis=0)
    for i in range(n_dims):
        if dot_products[i] < 0:
            new_volume = np.flip(new_volume, axis=i)
            aff_flo[:, i] = - aff_flo[:, i]
            aff_flo[:3, 3] = aff_flo[:3, 3] - aff_flo[:3, i] * (new_volume.shape[i] - 1)

    if return_aff:
        return new_volume, aff_flo
    else:
        return new_volume

def get_ras_axes(aff, n_dims=3):
    """This function finds the RAS axes corresponding to each dimension of a volume, based on its affine matrix.
    :param aff: affine matrix Can be a 2d numpy array of size n_dims*n_dims, n_dims+1*n_dims+1, or n_dims*n_dims+1.
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: two numpy 1d arrays of lengtn n_dims, one with the axes corresponding to RAS orientations,
    and one with their corresponding direction.
    """
    aff_inverted = np.linalg.inv(aff)
    img_ras_axes = np.argmax(np.absolute(aff_inverted[0:n_dims, 0:n_dims]), axis=0)
    for i in range(n_dims):
        if i not in img_ras_axes:
            unique, counts = np.unique(img_ras_axes, return_counts=True)
            incorrect_value = unique[np.argmax(counts)]
            img_ras_axes[np.where(img_ras_axes == incorrect_value)[0][-1]] = i

    return img_ras_axes


def unet(nb_features,
         input_shape,
         nb_levels,
         conv_size,
         nb_labels,
         name='unet',
         prefix=None,
         feat_mult=1,
         pool_size=2,
         padding='same',
         dilation_rate_mult=1,
         activation='elu',
         skip_n_concatenations=0,
         use_residuals=False,
         final_pred_activation='softmax',
         nb_conv_per_level=1,
         layer_nb_feats=None,
         conv_dropout=0,
         batch_norm=None,
         input_model=None):
    """
    Unet-style keras model with an overdose of parametrization. Copied with permission
    from github.com/adalca/neurite.
    """

    # naming
    model_name = name
    if prefix is None:
        prefix = model_name

    # volume size data
    ndims = len(input_shape) - 1
    if isinstance(pool_size, int):
        pool_size = (pool_size,) * ndims

    # get encoding model
    enc_model = conv_enc(nb_features,
                         input_shape,
                         nb_levels,
                         conv_size,
                         name=model_name,
                         prefix=prefix,
                         feat_mult=feat_mult,
                         pool_size=pool_size,
                         padding=padding,
                         dilation_rate_mult=dilation_rate_mult,
                         activation=activation,
                         use_residuals=use_residuals,
                         nb_conv_per_level=nb_conv_per_level,
                         layer_nb_feats=layer_nb_feats,
                         conv_dropout=conv_dropout,
                         batch_norm=batch_norm,
                         input_model=input_model)

    # get decoder
    # use_skip_connections=True makes it a u-net
    lnf = layer_nb_feats[(nb_levels * nb_conv_per_level):] if layer_nb_feats is not None else None
    dec_model = conv_dec(nb_features,
                         None,
                         nb_levels,
                         conv_size,
                         nb_labels,
                         name=model_name,
                         prefix=prefix,
                         feat_mult=feat_mult,
                         pool_size=pool_size,
                         use_skip_connections=True,
                         skip_n_concatenations=skip_n_concatenations,
                         padding=padding,
                         dilation_rate_mult=dilation_rate_mult,
                         activation=activation,
                         use_residuals=use_residuals,
                         final_pred_activation=final_pred_activation,
                         nb_conv_per_level=nb_conv_per_level,
                         batch_norm=batch_norm,
                         layer_nb_feats=lnf,
                         conv_dropout=conv_dropout,
                         input_model=enc_model)
    final_model = dec_model

    return final_model

def conv_enc(nb_features,
             input_shape,
             nb_levels,
             conv_size,
             name=None,
             prefix=None,
             feat_mult=1,
             pool_size=2,
             dilation_rate_mult=1,
             padding='same',
             activation='elu',
             layer_nb_feats=None,
             use_residuals=False,
             nb_conv_per_level=2,
             conv_dropout=0,
             batch_norm=None,
             input_model=None):
    """
    Fully Convolutional Encoder. Copied with permission from github.com/adalca/neurite.
    """

    # naming
    model_name = name
    if prefix is None:
        prefix = model_name

    # first layer: input
    name = '%s_input' % prefix
    if input_model is None:
        input_tensor = keras.layers.Input(shape=input_shape, name=name)
        last_tensor = input_tensor
    else:
        input_tensor = input_model.inputs
        last_tensor = input_model.outputs
        if isinstance(last_tensor, list):
            last_tensor = last_tensor[0]

    # volume size data
    ndims = len(input_shape) - 1
    if isinstance(pool_size, int):
        pool_size = (pool_size,) * ndims

    # prepare layers
    convL = getattr(keras.layers, 'Conv%dD' % ndims)
    conv_kwargs = {'padding': padding, 'activation': activation, 'data_format': 'channels_last'}
    maxpool = getattr(keras.layers, 'MaxPooling%dD' % ndims)

    # down arm:
    # add nb_levels of conv + ReLu + conv + ReLu. Pool after each of first nb_levels - 1 layers
    lfidx = 0  # level feature index
    for level in range(nb_levels):
        lvl_first_tensor = last_tensor
        nb_lvl_feats = np.round(nb_features * feat_mult ** level).astype(int)
        conv_kwargs['dilation_rate'] = dilation_rate_mult ** level

        for conv in range(nb_conv_per_level):  # does several conv per level, max pooling applied at the end
            if layer_nb_feats is not None:  # None or List of all the feature numbers
                nb_lvl_feats = layer_nb_feats[lfidx]
                lfidx += 1

            name = '%s_conv_downarm_%d_%d' % (prefix, level, conv)
            if conv < (nb_conv_per_level - 1) or (not use_residuals):
                last_tensor = convL(nb_lvl_feats, conv_size, **conv_kwargs, name=name)(last_tensor)
            else:  # no activation
                last_tensor = convL(nb_lvl_feats, conv_size, padding=padding, name=name)(last_tensor)

            if conv_dropout > 0:
                # conv dropout along feature space only
                name = '%s_dropout_downarm_%d_%d' % (prefix, level, conv)
                noise_shape = [None, *[1] * ndims, nb_lvl_feats]
                last_tensor = keras.layers.Dropout(conv_dropout, noise_shape=noise_shape, name=name)(last_tensor)

        if use_residuals:
            convarm_layer = last_tensor

            # the "add" layer is the original input
            # However, it may not have the right number of features to be added
            nb_feats_in = lvl_first_tensor.get_shape()[-1]
            nb_feats_out = convarm_layer.get_shape()[-1]
            add_layer = lvl_first_tensor
            if nb_feats_in > 1 and nb_feats_out > 1 and (nb_feats_in != nb_feats_out):
                name = '%s_expand_down_merge_%d' % (prefix, level)
                last_tensor = convL(nb_lvl_feats, conv_size, **conv_kwargs, name=name)(lvl_first_tensor)
                add_layer = last_tensor

                if conv_dropout > 0:
                    name = '%s_dropout_down_merge_%d_%d' % (prefix, level, conv)
                    noise_shape = [None, *[1] * ndims, nb_lvl_feats]

            name = '%s_res_down_merge_%d' % (prefix, level)
            last_tensor = keras.layers.add([add_layer, convarm_layer], name=name)

            name = '%s_res_down_merge_act_%d' % (prefix, level)
            last_tensor = keras.layers.Activation(activation, name=name)(last_tensor)

        if batch_norm is not None:
            name = '%s_bn_down_%d' % (prefix, level)
            last_tensor = keras.layers.BatchNormalization(axis=batch_norm, name=name)(last_tensor)

        # max pool if we're not at the last level
        if level < (nb_levels - 1):
            name = '%s_maxpool_%d' % (prefix, level)
            last_tensor = maxpool(pool_size=pool_size, name=name, padding=padding)(last_tensor)

    # create the model and return
    model = keras.Model(inputs=input_tensor, outputs=[last_tensor], name=model_name)
    return model


def conv_dec(nb_features,
             input_shape,
             nb_levels,
             conv_size,
             nb_labels,
             name=None,
             prefix=None,
             feat_mult=1,
             pool_size=2,
             use_skip_connections=False,
             skip_n_concatenations=0,
             padding='same',
             dilation_rate_mult=1,
             activation='elu',
             use_residuals=False,
             final_pred_activation='softmax',
             nb_conv_per_level=2,
             layer_nb_feats=None,
             batch_norm=None,
             conv_dropout=0,
             input_model=None):
    """
    Fully Convolutional Decoder. Copied with permission from github.com/adalca/neurite.

    Parameters:
        ...
        use_skip_connections (bool): if true, turns an Enc-Dec to a U-Net.
            If true, input_tensor and tensors are required.
            It assumes a particular naming of layers. conv_enc...
    """

    # naming
    model_name = name
    if prefix is None:
        prefix = model_name

    # if using skip connections, make sure need to use them.
    if use_skip_connections:
        assert input_model is not None, "is using skip connections, tensors dictionary is required"

    # first layer: input
    input_name = '%s_input' % prefix
    if input_model is None:
        input_tensor = keras.layers.Input(shape=input_shape, name=input_name)
        last_tensor = input_tensor
    else:
        input_tensor = input_model.input
        last_tensor = input_model.output
        input_shape = last_tensor.shape.as_list()[1:]

    # vol size info
    ndims = len(input_shape) - 1
    if isinstance(pool_size, int):
        if ndims > 1:
            pool_size = (pool_size,) * ndims

    # prepare layers
    convL = getattr(keras.layers, 'Conv%dD' % ndims)
    conv_kwargs = {'padding': padding, 'activation': activation}
    upsample = getattr(keras.layers, 'UpSampling%dD' % ndims)

    # up arm:
    # nb_levels - 1 layers of Deconvolution3D
    #    (approx via up + conv + ReLu) + merge + conv + ReLu + conv + ReLu
    lfidx = 0
    for level in range(nb_levels - 1):
        nb_lvl_feats = np.round(nb_features * feat_mult ** (nb_levels - 2 - level)).astype(int)
        conv_kwargs['dilation_rate'] = dilation_rate_mult ** (nb_levels - 2 - level)

        # upsample matching the max pooling layers size
        name = '%s_up_%d' % (prefix, nb_levels + level)
        last_tensor = upsample(size=pool_size, name=name)(last_tensor)
        up_tensor = last_tensor

        # merge layers combining previous layer
        if use_skip_connections & (level < (nb_levels - skip_n_concatenations - 1)):
            conv_name = '%s_conv_downarm_%d_%d' % (prefix, nb_levels - 2 - level, nb_conv_per_level - 1)
            cat_tensor = input_model.get_layer(conv_name).output
            name = '%s_merge_%d' % (prefix, nb_levels + level)
            last_tensor = keras.layers.concatenate([cat_tensor, last_tensor], axis=ndims + 1, name=name)

        # convolution layers
        for conv in range(nb_conv_per_level):
            if layer_nb_feats is not None:
                nb_lvl_feats = layer_nb_feats[lfidx]
                lfidx += 1

            name = '%s_conv_uparm_%d_%d' % (prefix, nb_levels + level, conv)
            if conv < (nb_conv_per_level - 1) or (not use_residuals):
                last_tensor = convL(nb_lvl_feats, conv_size, **conv_kwargs, name=name)(last_tensor)
            else:
                last_tensor = convL(nb_lvl_feats, conv_size, padding=padding, name=name)(last_tensor)

            if conv_dropout > 0:
                name = '%s_dropout_uparm_%d_%d' % (prefix, level, conv)
                noise_shape = [None, *[1] * ndims, nb_lvl_feats]
                last_tensor = keras.layers.Dropout(conv_dropout, noise_shape=noise_shape, name=name)(last_tensor)

        # residual block
        if use_residuals:

            # the "add" layer is the original input
            # However, it may not have the right number of features to be added
            add_layer = up_tensor
            nb_feats_in = add_layer.get_shape()[-1]
            nb_feats_out = last_tensor.get_shape()[-1]
            if nb_feats_in > 1 and nb_feats_out > 1 and (nb_feats_in != nb_feats_out):
                name = '%s_expand_up_merge_%d' % (prefix, level)
                add_layer = convL(nb_lvl_feats, conv_size, **conv_kwargs, name=name)(add_layer)

                if conv_dropout > 0:
                    name = '%s_dropout_up_merge_%d_%d' % (prefix, level, conv)
                    noise_shape = [None, *[1] * ndims, nb_lvl_feats]
                    last_tensor = keras.layers.Dropout(conv_dropout, noise_shape=noise_shape, name=name)(last_tensor)

            name = '%s_res_up_merge_%d' % (prefix, level)
            last_tensor = keras.layers.add([last_tensor, add_layer], name=name)

            name = '%s_res_up_merge_act_%d' % (prefix, level)
            last_tensor = keras.layers.Activation(activation, name=name)(last_tensor)

        if batch_norm is not None:
            name = '%s_bn_up_%d' % (prefix, level)
            last_tensor = keras.layers.BatchNormalization(axis=batch_norm, name=name)(last_tensor)

    # Compute likelyhood prediction (no activation yet)
    name = '%s_likelihood' % prefix
    last_tensor = convL(nb_labels, 1, activation=None, name=name)(last_tensor)
    like_tensor = last_tensor

    # output prediction layer
    # we use a softmax to compute P(L_x|I) where x is each location
    if final_pred_activation == 'softmax':
        name = '%s_prediction' % prefix
        softmax_lambda_fcn = lambda x: keras.activations.softmax(x, axis=ndims + 1)
        pred_tensor = keras.layers.Lambda(softmax_lambda_fcn, name=name)(last_tensor)

    # otherwise create a layer that does nothing.
    else:
        name = '%s_prediction' % prefix
        pred_tensor = keras.layers.Activation('linear', name=name)(like_tensor)

    # create the model and retun
    model = keras.Model(inputs=input_tensor, outputs=pred_tensor, name=model_name)
    return model


def save_volume(volume, aff, header, path, res=None, dtype=None, n_dims=3):
    """
    Save a volume.
    :param volume: volume to save
    :param aff: affine matrix of the volume to save. If aff is None, the volume is saved with an identity affine matrix.
    aff can also be set to 'FS', in which case the volume is saved with the affine matrix of FreeSurfer outputs.
    :param header: header of the volume to save. If None, the volume is saved with a blank header.
    :param path: path where to save the volume.
    :param res: (optional) update the resolution in the header before saving the volume.
    :param dtype: (optional) numpy dtype for the saved volume.
    :param n_dims: (optional) number of dimensions, to avoid confusion in multi-channel case. Default is None, where
    n_dims is automatically inferred.
    """

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


# execute script
if __name__ == '__main__':
    main()
