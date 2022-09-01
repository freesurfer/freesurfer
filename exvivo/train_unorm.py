#!/usr/bin/env python
# coding: utf-8

# # synthseg quick start (tf 2.0+)
# 
# `PYTHONPATH` should include a path to several libraries, including:  
# `neurite`, `voxelmorph`, `neurite-sandbox`, `voxelmorph-sandbox`, `pystrum`
# 
# 
# 


import socket, os
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers as KL
from tqdm import tqdm
#import glob, copy
#import nibabel as nib
from pathlib import Path
import tensorflow_probability as tfp
from tensorflow.keras import backend as K

import freesurfer as fs
#from freesurfer import deeplearn as fsd

import neurite as ne
import neurite_sandbox as nes
#import voxelmorph as vxm
#import voxelmorph_sandbox as vxms
from neurite_sandbox.tf.utils.utils import plot_fit_callback as pfc
#from tensorflow.compat.v1 import ConfigProto
#from tensorflow.compat.v1 import InteractiveSession


sigma = 1
dofit = False
dotrans = True
which_loss = 'mse'
which_opt = 'adam'
lsigma = False
last_nfilters=0
use_exp = True
fit_bias = True
only_wm = True
save_model = True
power = 2
bias_wt = .01

hemi = 'lh'
"""hemi = 'rh'"""
def erase_labels(seg, img, erase_list):
    img = img.copy()
    img[np.isin(seg, erase_list)] = 0
    return img

print(f'host name {socket.gethostname()}')
print(f'dofit is {dofit}, dotrans is {dotrans}, which_loss is {which_loss}, which_opt {which_opt}, lsigma is {lsigma}, exp is {use_exp}, hemi {hemi}, only_wm {only_wm}, sigma {sigma}, bias_wt {bias_wt}, fit_bias {fit_bias}, save_model {save_model}')


gpu_id = 0
os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

if 'inited' not in locals():
    #tf.config.threading.set_intra_op_parallelism_threads(6)
    #tf.config.threading.set_inter_op_parallelism_threads(6)
    inited = 1

#fsd.configure(gpu=gpu_id)
ne.utils.setup_device(gpuid=gpu_id)
#config = ConfigProto()
#config.gpu_options.allow_growth = True
#session = InteractiveSession(config=config)

policy = tf.keras.mixed_precision.experimental.Policy('mixed_float16')
os.environ['TF_XLA_FLAGS'] = '--tf_xla_auto_jit=2 --tf_xla_cpu_global_jit'
tf.compat.v1.enable_eager_execution()


odir = '/autofs/cluster/vxmdata1/FS_Slim/proc/oasis/dataset'
#lut = fs.lookups.LookupTable.read(os.path.join(odir, 'seg4_labels.txt'))
lut = fs.lookups.LookupTable.read(os.path.join(odir, 'seg35_labels.txt'))

nlabels = len(lut)

# get subjects
subjects = [f for f in Path(odir).iterdir() if 'OASIS' in str(f)]

# read in label maps
target_shape = (256,)*3
seg_files = [f/f'seg{nlabels-2}.nii.gz' for f in subjects]
ntest = 10
crop = -1 if dofit else ntest
seg_files = seg_files[:crop]
lh_list = []
rh_list = []
bs_list = []
bs_names = [
    'Brain-Stem',
    '4th-Ventricle',
    'Left-Cerebellum-White-Matter',
    'Left-Cerebellum-Cortex',
    'Right-Cerebellum-White-Matter',
    'Right-Cerebellum-Cortex'
]
if 'seg_list' not in locals() and 'seg_list' not in globals():
    for ind in lut.keys():
        if lut[ind].name in bs_names:
            bs_list.append(ind)
        if lut[ind].name.startswith('Left'):
            lh_list.append(ind)
        elif lut[ind].name.startswith('Right'):
            rh_list.append(ind)

    seg_list = [fs.Volume.read(str(f)).fit_to_shape(target_shape, center='bbox') for f in tqdm(seg_files)]
    if hemi == 'lh':
        erase_list = rh_list
    else:
        erase_list = lh_list

    for seg in tqdm(seg_list):
        seg.data = erase_labels(seg.data, seg.data, erase_list)

    segs = np.array([mri.data[...,np.newaxis] for mri in seg_list])
    bag_label = np.array(segs).max()+1
    if 0:
        lab_to_ind0, ind_to_lab0 = ne.py.utils.rebase_lab(np.array(segs))
        bag_label = len(lab_to_ind0)
        lab_to_ind = np.zeros(len(lab_to_ind0)+1,).astype(int)
        ind_to_lab = np.zeros(len(ind_to_lab0)+1,).astype(int)
        lab_to_ind[:-1] = lab_to_ind0
        ind_to_lab[:-1] = ind_to_lab0
        lab_to_ind[bag_label] = len(ind_to_lab0)
        ind_to_lab[-1] = bag_label

vol_shape = list(seg_list[0].shape)

batch_size = 1
border = 5
nfeatures = 80 // batch_size
min_shape = np.min(vol_shape)

nb_conv_per_level = 3
nfeatures = 96 // batch_size

for max_nlevels in range(1,np.floor(np.log2(min_shape)).astype(int)):
    if 2**max_nlevels > min_shape:
        break

nlevels = 4
ndec = 4
ndown = 0

# pass seg to model but only so it can be used as a mask in output
seg_input  = KL.Input(vol_shape + [1])
img_input  = KL.Input(vol_shape + [1])
if ndown > 0:
    zooms = [2**ndown]*3
    down_tensor = ne.layers.Resize(list(1./np.array(zooms)), name=f'down_tensor')(img_input)
    down_shape = list(np.array(target_shape) // 2**ndown)
    unet = ne.models.unet(nfeatures, down_shape + [1], nlevels, 3, 1, name='inorm', nb_conv_per_level=nb_conv_per_level, final_pred_activation='linear')
    unet_output = unet(down_tensor)
    up_tensor = ne.layers.Resize(zooms, name=f'up_tensor')(unet_output)
else:
    nf = nb_conv_per_level * [[16]] + nlevels * [nb_conv_per_level * [nfeatures]]
    unet = ne.models.unet(nf, vol_shape + [1], None, 3, 1, name='inorm', nb_conv_per_level=None, final_pred_activation='linear', feat_mult=None)
    up_tensor = unet(img_input)

if last_nfilters > 0:
    convL = getattr(KL, 'Conv%dD' % 3)
    smoothed = convL(last_nfilters, kernel_size=3, padding='same', name=f'smoothed{last_nfilters}')(up_tensor)
    if lsigma:
        smoothed = convL(1, kernel_size=3, padding='same', name='conv_out')(smoothed)
    else:
        smoothed = convL(1, kernel_size=3, padding='same', name='conv_out')(smoothed)
        smoothed = nes.layers.GaussianBlur(sigma=sigma, name='conv_out_smoothed')(smoothed)
else:
    if sigma > 0:
        smoothed = nes.layers.GaussianBlur(sigma=sigma, name='conv_out_smoothed')(up_tensor)
    else:
        smoothed = up_tensor
    
if use_exp:
    f = lambda x : tf.exp(x)
else:
    f = lambda x: tf.clip_by_value(1. + x, clip_value_min=0, clip_value_max=4)

if which_loss == 'mse':  # don't have to worry about making the field mean 1
    bias_field = KL.Lambda(f, name='bias_field')(smoothed)
else:
    unorm_field = KL.Lambda(f, name='unorm_field')(smoothed)
    bias_mean = KL.Lambda(lambda x: K.mean(x, keepdims=True)-1)(unorm_field)
    bias_field = KL.Subtract(name='bias_field')([unorm_field, bias_mean])

img_out = KL.Multiply(name='img_corrected')([img_input, bias_field])
smodel = tf.keras.Model(img_input, img_out)

smodel_out = smodel.outputs[0]
nrounds = 1
if nrounds > 1:
    for round in range(1, nrounds):
        smodel_out = smodel(smodel_out)

    # model = tf.keras.Model(inp_tensor, model_out)

model_out = KL.Concatenate(name='out_concat')([smodel_out, seg_input])
outputs = [model_out]
if fit_bias:
    bias_stack = KL.Concatenate(name='bias_stack')([bias_field, seg_input])
    outputs = outputs + [bias_stack]
model = tf.keras.Model([img_input, seg_input], outputs)

print(f'encoder with {nlevels} levels, {nrounds} rounds')

model.summary(line_length=120)

# labels to image variant
mean_max = [110] * (nlabels+1)
mean_max[0] = 40
bias_res = 1
bias_res_list = (12, 16, 32, 64, 96) if bias_res == 2 else (16, 32, 64)
bias_scale = .5
#nes.models.labels_to_image(..., bias_exponential=False)

clip_func = lambda x: tf.clip_by_value(1. + x, clip_value_min=0, clip_value_max=bias_max * 10)
gen_arg_1 = dict(
    in_shape=vol_shape,
    std_max=15,
    zero_background=.8,
    out_shape=vol_shape,
    labels_in=range(nlabels+1),   # +1 is for bag
    one_hot=False,
    warp_max=.5*6,
    warp_res=64,
    #bias_res=bias_res * 16,   # was 16, then 4*16
    # bias_std=bias_scale*2*.5,   # was .5 then 2*.5
    bias_res=bias_res_list,
    bias_max=bias_scale,
    max_shift=20*dotrans,
    max_rotate=15*dotrans,
    max_shear=.2*dotrans,
    max_scale=.3*dotrans,
    dc_offset=40,
    # gamma_std=.25,
    clip_max=2800,
    normalize=False,
    mean_max=mean_max,
    return_mean=True,
    return_bias=True
)
gen_arg_2 = gen_arg_1.copy()
gen_arg_2['bias_std'] = 0
tf.random.set_seed(1)
# noise
seeds = dict(mean=1, std=2, warp=3, blur=4, bias=5, gamma=6, shift=7, rotate=8, scale=9, shear=10, dc_offset=11, background=12)
use_exp_in_gen = False
if use_exp:
    f = tf.exp
else:
    f = lambda x: tf.clip_by_value(1. + x, clip_value_min=0, clip_value_max=4)
#gen_model_1 = nes.models.labels_to_image(**gen_arg_1, id=1, seeds=seeds, bias_func=f)
std_func=nes.utils.generative.sample_identical
gen_model_1 = nes.models.labels_to_image(**gen_arg_1, id=1, seeds=seeds,
                                         std_func=std_func, bias_func=f)
# gen_model_2 = nes.models.labels_to_image(**gen_arg_2, id=2, seeds=seeds)

        
def bias_func(x, bias_prob=0.5):
    # Bias-conditioning function of choice.
    x = tf.exp(x)

    # Randomized application.
    if bias_prob < 1:
        bias_bit = tf.less(tf.random.uniform(shape=[]), bias_prob)
        bias_not = tf.logical_not(bias_bit)
        x = x * tf.cast(bias_bit, x.dtype) + tf.cast(bias_not, x.dtype)

    return x


@tf.autograph.experimental.do_not_convert
def generator(segs, lh_list, rh_list, bs_list, 
              lh_keep_pval=.5, rh_keep_pval=.5, bs_keep_pval=1, 
              insert_bag_pval=.8, batch_size=16,
              max_bag_dilations=25, use_rand=True, 
              seg_dilations=0, seg_closes=0, lab_to_ind=None,
              return_bias=False):

    bag_label = np.array(segs).max()+1

    if not use_rand:
        idx0 = 0
    while 1:
        if not use_rand:
            idx = np.array([np.mod(ind, segs.shape[0]) for ind in range(idx0, idx0+batch_size)])
            idx0 = np.mod(idx0+batch_size, segs.shape[0])
        else:
            idx = np.random.randint(0, segs.shape[0], batch_size)

        seg_in = segs[idx, ...]
        erase_list = []
        if use_rand:
            if np.random.randint(2):  # check lh first
                if np.random.rand() > lh_keep_pval:
                    erase_list += lh_list
                elif np.random.rand() > rh_keep_pval:
                    erase_list += rh_list
            else:    # check rh first to keep it unbiased
                if np.random.rand() > rh_keep_pval:
                    erase_list += rh_list
                elif np.random.rand() > lh_keep_pval:
                    erase_list += lh_list

            if np.random.rand() > bs_keep_pval:
                erase_list += bs_list
        else:  # disable randomness for testing/debugging
            if (np.mod(idx0,3) and lh_keep_pval < 1) or (lh_keep_pval == 0):
                erase_list += lh_list
            elif  (np.mod(idx0, 5) and rh_keep_pval < 1) or (rh_keep_pval == 0):
                erase_list += rh_list
            if (np.mod(idx0,7) and bs_keep_pval < 1) or (bs_keep_pval == 0):
                erase_list += bs_list

        seg_in = erase_labels(seg_in, seg_in, erase_list)

        brain_mask = tf.cast(tf.where(seg_in > 0, 1, 0), tf.float32)
        if np.random.rand() < insert_bag_pval:  # put bag in sometimes
            bag_dilations = np.random.randint(1, max_bag_dilations)

            crop_lims = [26 /27.0, 27 / 27.0]
            crop_lims = [0, .5]
            dil_fn1 = lambda x: nes.utils.utils.morphology_3d(x, 1, 1, 
                                                             operation='dilate', 
                                                             eight_connectivity=False)
            seg_dil = tf.map_fn(dil_fn1, brain_mask, fn_output_signature=tf.bool)
            dil_fn2 = lambda x: nes.utils.utils.morphology_3d(x, 1, bag_dilations, 
                                                             operation='dilate', 
                                                             rand_crop=crop_lims,
                                                             eight_connectivity=False)
            for i in range(3):
                seg_dil = tf.cast(seg_dil, dtype=brain_mask.dtype)
                seg_dil = tf.map_fn(dil_fn2, seg_dil, fn_output_signature=tf.bool)
            seg_bag = tf.cast(seg_dil, brain_mask.dtype) - brain_mask
            seg_in += tf.cast(seg_bag, seg_in.dtype) * (bag_label)  # new label

        img, seg, img_no_bias,bias = gen_model_1.predict(seg_in)

        if lab_to_ind is not None:
            seg_old = seg.copy()
            seg = lab_to_ind[seg]
            # gdb.set_trace()

        if not use_rand and idx0 == 1:
            # gdb.set_trace()
            idx0 *= 1

        #img_no_bias, seg, mean_img2 = gen_model_2.predict(seg_in)
        low_mask = tf.cast(seg > 0, img.dtype)
        hi_mask = tf.cast(seg < bag_label, img.dtype)
        brain_mask = low_mask * hi_mask
        if seg_closes > 0:
            dil_fn = lambda x: nes.utils.utils.morphology_3d(x, 1, seg_closes, 
                                                             operation='close', 
                                                             eight_connectivity=True)
            brain_mask = tf.map_fn(dil_fn, brain_mask, fn_output_signature=tf.bool).numpy().astype('float32')

        brain_vol = brain_mask * img
        brain_vol_no_bias = img_no_bias * brain_mask

        # map the median in the unnormalized image to 0.5 in the normed one
        #median = tf.cast(tfp.stats.percentile(brain_vol[brain_vol>0], 50), img.dtype)
        #norm_val = median / .5
        # tu = tf.unique_with_counts(tf.cast(brain_vol[brain_mask>0], tf.int16))
        norm_val = tf.cast(tfp.stats.percentile(brain_vol[brain_mask>0], 99), img.dtype)
        # norm_val = tf.cast(tu[0][0], tf.float32) / .5  # map most common brain value to .5
        img = tf.clip_by_value(tf.math.divide_no_nan(img, norm_val), 0, 3)
        img_no_bias = tf.clip_by_value(tf.math.divide_no_nan(img_no_bias, norm_val), 0, 3)

        # now make the image means the same
        brain_vol = brain_mask * img
        brain_vol_no_bias = brain_mask * img_no_bias
        img_mean = tf.reduce_mean(brain_vol[brain_mask > 0])
        img_no_bias_mean = tf.reduce_mean(brain_vol_no_bias[brain_mask > 0])

        ratio = tf.math.divide_no_nan(img_no_bias_mean, img_mean)
        img *= ratio
        # img *= brain_mask
        img_no_bias *= brain_mask
        # seg *= tf.cast(brain_mask, dtype=seg.dtype)
        #inputs = np.concatenate([img.numpy(), seg], axis=-1)
        outputs = img_no_bias.numpy()
        if return_bias:
            outputs = [outputs, 1./(bias+1e-9)]
        yield [img.numpy(), seg], outputs



bs_keep_pval = .8
insert_bag_pval = .8
if hemi == 'lh':
    lh_keep_pval = 1
    rh_keep_pval = 0
else:
    lh_keep_pval = 0
    rh_keep_pval = 1

gen = generator(segs, lh_list, rh_list, bs_list, batch_size=batch_size,
                lh_keep_pval=lh_keep_pval, rh_keep_pval=rh_keep_pval, 
                bs_keep_pval=bs_keep_pval, insert_bag_pval=insert_bag_pval,
                seg_closes=6, lab_to_ind=None, return_bias=fit_bias)


import sys
class label_coefficient_of_variation:
    def __init__(self, label_weights=None, power=2):
        self.label_weights = tf.convert_to_tensor(label_weights, tf.float32)

    ''' 
    compute the mean and variance for each label separately, then average over those so that
    the classes are equally weighted regardless of frequency

    returns the coefficient of variation of each class in a (1,nlabels-1) tensor
    the bg label is ignored (which is why it is nlabels-1)

    y_pred should have a ground-truth label map in frame 1 
    '''
    def loss(self, y_true, y_pred):
        label_map = tf.keras.utils.to_categorical(y_pred[...,1:2])
        pred_intensity_map = y_pred[...,0:1]
        
        inshape = tf.shape(y_true)
        ndim = len(inshape) - 2
        nlabels = tf.shape(label_map)[-1]
        pred_intensity_vec = tf.reshape(pred_intensity_map, [1, -1])    # 1 x nvoxels
        
        label_vec = tf.reshape(label_map, [-1, nlabels])  # nvoxels x nlabels
        label_sum_sq = tf.linalg.matmul(tf.math.square(pred_intensity_vec), label_vec)  # 1 x nlabels
        pred_label_sum = tf.linalg.matmul(pred_intensity_vec, label_vec)  # 1 x nlabels
        label_nvox = tf.reduce_sum(label_vec, axis=0)        # number of vox in each label
        pred_label_mean = tf.math.divide_no_nan(pred_label_sum, label_nvox)
        pred_label_mean_sq = tf.square(pred_label_mean)
        label_var = tf.math.divide_no_nan(label_sum_sq, label_nvox) - pred_label_mean_sq
        if power > 2:
            pred_label_mean_sq = tf.pow(pred_label_mean, power)
            label_var = tf.pow(label_var, power-2)
            
        
        # [:,1:] ignores background voxels (where label is 0). Make sure 0/0 goes to 1 not 0
        eps = sys.float_info.epsilon
        if self.label_weights is not None:   # keep bg label in and allow user to decide on weighting
            # cut off bag label for vols that don't have it
            label_weights = self.label_weights[:nlabels] 
            numer = label_var + eps
            denom = pred_label_mean_sq + eps
            scale = np.prod(np.array(label_weights.get_shape().as_list()[1:])) / tf.reduce_sum(label_weights)
            coef_var = scale * tf.math.divide_no_nan(label_weights*numer, denom)
        else:
            numer = label_var[:,1:] + eps
            denom = pred_label_mean_sq[:,1:] + eps
            coef_var = tf.math.divide_no_nan(numer, denom)

        return coef_var


class masked_mse_loss:
    def __init__(self, keep_labels=None, debug=False):
        self.debug = debug
        if keep_labels is not None and type(keep_labels) is not list:
            self.keep_labels = [keep_labels]
        else:
            self.keep_labels = keep_labels

    def loss(self, y_true, y_pred):
        eps = sys.float_info.epsilon
        if self.keep_labels is not None:
            mask = tf.zeros(y_true.shape, dtype=tf.bool)
            for keep_label  in self.keep_labels:
                mask = tf.logical_or(mask, y_pred[...,1:2] == keep_label)

            mask = tf.cast(mask, y_pred.dtype)
        else:   # if none specified, just use all non-zero labels
            mask = tf.cast(y_pred[...,1:2] > 0, y_pred.dtype)
        pred_image = tf.math.log(y_pred[...,0:1]+eps)
        true_image = tf.math.log(y_true + eps)
        dif = tf.math.squared_difference(true_image, pred_image)
        if self.debug:
            fv = fs.Freeview(swap_batch_dim=True)
            fv.vol(true_image, name='true', opts=':linked=1:visible=0')
            fv.vol(pred_image, name='pred', opts=':linked=1')
            fv.vol(dif*mask, name='dif', opts=':linked=1:colormap=heat:locked=1')
            fv.show()
            
        lval = 10. * tf.math.divide_no_nan(tf.reduce_sum(dif * mask), tf.reduce_sum(mask))
        return lval


class weighted_loss:
    def __init__(self, var_wt=1, mse_wt=1):
        self.var_wt = var_wt
        self.mse_wt = mse_wt

    def loss(self, y_true, y_pred):
        var_loss = label_variance(y_true, y_pred) 
        mse_loss = masked_mse_loss(y_true, y_pred)
        return var_loss * self.var_wt + mse_loss * self.mse_wt


scale = 1
epsilon = 1e-6
lr = 1e-5
warm = 0

if warm == 1:
    #warm_restart_epoch=[5, 10, 15, 20]
    #warm_restart_lr=[1e-3, 2e-3, 5e-3, 1e-2]
    #warm_restart_epoch=[5]
    #warm_restart_lr=[1e-3]
    warm_restart_epoch=[10, 20, 50]
    warm_restart_lr=[2*lr, 5*lr, 10*lr]
elif warm == 2:
    warm_restart_epoch=[5, 10, 15, 20, 150]
    warm_restart_lr=[1e-3, 2e-3, 5e-3, 1e-2, 5e-2]
    warm_restart_epoch=[5, 10, 15, 20, 150]
    warm_restart_lr=[1e-3, 2e-3, 5e-3, 1e-2, 1e-2]
    warm_restart_epoch=[50]
    warm_restart_lr=[lr]
elif warm == 3:
    warm_restart_epoch=[5, 10, 20]
    warm_restart_lr=[2*lr, 5*lr, 10*lr]
elif warm == 4:
    warm_restart_epoch=[5, 10, 15]
    warm_restart_lr=[1e-3, 2e-3, 5e-3]
elif warm == 5:
    warm_restart_epoch=[5, 10, 15]
    warm_restart_lr=[1e-3, 2e-3, 3e-3]
else:
    warm_restart_epoch = None
    warm_restart_lr = None

debug = False

name = f'{hemi}.exvivo.unet.norm.nf{nfeatures}.nl{nlevels}.trans.{dotrans}.conv_per_level.{nb_conv_per_level}.lr.{lr:.1e}.eps.{epsilon:.1e}.warm{warm}.lsigma.{lsigma}.sigma.{sigma}.ndown.{ndown}.loss.{which_loss}.bias.{bias_scale}.bres.{bias_res}.exp.{use_exp}.nrounds.{nrounds}.nlast.{last_nfilters}.only_wm.{only_wm}.fit_bias.{fit_bias}'
if bias_wt > 0:
    name = f'{hemi}.exvivo.unet.norm.nf{nfeatures}.nl{nlevels}.trans.{dotrans}.conv_per_level.{nb_conv_per_level}.lr.{lr:.1e}.eps.{epsilon:.1e}.warm{warm}.lsigma.{lsigma}.sigma.{sigma}.ndown.{ndown}.loss.{which_loss}.bias.{bias_scale}.bres.{bias_res}.exp.{use_exp}.nrounds.{nrounds}.nlast.{last_nfilters}.only_wm.{only_wm}.fit_bias.{bias_wt}'
if debug:
    name += '.debug'

if which_loss == 'mse':
    search_str = 'Left' if hemi == 'lh' else 'Right'
    search_str += '-Cerebral-'
    if only_wm:
        search_str += 'White-Matter'
    
    keep_labels = lut.search(search_str)
    loss = masked_mse_loss(keep_labels=keep_labels, debug=False).loss
    loss_bias = masked_mse_loss(debug=False).loss
elif which_loss == 'combo':
    loss = weighted_loss(var_wt=.01).loss
elif which_loss == 'cov':
    if only_wm:
        label_weights = np.zeros((bag_label+1,))
        label_weights[1] = 1   # cerebral wm 
        #label_weights[5] = 1   # cerebellar wm
        #label_weights[13] = 1  # brainstem
    else:
        label_weights = np.ones((bag_label+1,))
        label_weights[0] = 0
        label_weights[lab_to_ind[bag_label]] = 0

    loss = label_coefficient_of_variation(label_weights, power=power).loss

if debug:
    losses = [nes.losses.debug_loss(loss, thresh=50, burn_in_epochs=10, input_stacked=False).loss]
else:
    losses = [loss]

if fit_bias:
    losses += [loss_bias]

loss_weights = [1]*len(model.outputs)
if fit_bias:
    loss_weights[1] = bias_wt
    loss_weights[0] = 1-bias_wt

opt = tf.keras.optimizers.Adam(learning_rate=lr, epsilon=epsilon) 
nes.utils.check_and_compile(model, gen,
                            optimizer=opt,
                            loss=losses,
                            loss_weights=loss_weights,
                            extra_outputs=[1, 1],
                            check_layers=False, run_eagerly=True)

initial_epoch = 0
if initial_epoch > 0:
    write_mode = 'a'
    fname = name + '.h5'
    print(f'loading weights from {fname}')
    model.load_weights(fname)

    if 0:
        lr_epoch = lr
        lr_scale = .5
        warm_step = 100
        warm_restart_epoch = []
        warm_restart_lr = []
        for epoch in range(initial_epoch+1, (initial_epoch+1)+11*warm_step, warm_step):
            lr_epoch += lr_scale*lr
            warm_restart_epoch += [epoch]
            warm_restart_lr += [lr_epoch]

        print(f'epochs {warm_restart_epoch}, lr {warm_restart_lr}')
else:
    write_mode = 'w'

write_cb = nes.callbacks.WriteHist(name+f'.txt', mode=write_mode)
lr_cb = nes.tf.callbacks.ReduceLRWithModelCheckpointAndRecovery(
    name+'.h5', 
    monitor='loss',
    verbose=1, 
    cooldown=60, 
    factor=.8, 
    thresh=None if initial_epoch == 0  else 2,
    patience=600, 
    thresh_inc_factor=6,  # 1.5,
    save_weights_only=True, 
    min_lr=1e-7, 
    burn_in=20,
    recovery_decrease_factor=1,
    warm_restart_epoch=warm_restart_epoch,
    nloss=5,
    warm_restart_lr=warm_restart_lr)

callbacks = [lr_cb, write_cb]

if dofit:
    print(f'saving fit results to {write_cb.fname} and model to {lr_cb.fname}')
    hist = model.fit(gen, epochs=10000, steps_per_epoch=50, callbacks=callbacks, 
                     initial_epoch=initial_epoch)
#                     initial_epoch=initial_epoch, use_multiprocessing=True, workers=6)
else:
    if 1:
        model.load_weights(lr_cb.fname)

test_gen = generator(segs, lh_list, rh_list, bs_list, batch_size=batch_size,
                     lh_keep_pval=lh_keep_pval, rh_keep_pval=rh_keep_pval,
                     bs_keep_pval=bs_keep_pval, insert_bag_pval=insert_bag_pval,
                     use_rand=False, seg_closes=6, return_bias=True, max_bag_dilations=25)

#test_gen = generator(segs, lh_list, rh_list, bs_list, batch_size=1,
#                     lh_keep_pval=.5, rh_keep_pval=.5, bs_keep_pval=.75, use_rand=False)
#                     
inb, outb = next(gen)

print('model eval:')
model.evaluate(inb, outb)

pl = []
pil = []
il = []
ol = []
sl = []
bl = []
bt = []
if fit_bias:
    bmodel = model
else:
    bias_tensor = model.get_layer('bias_field').output
    bmodel = tf.keras.Model(model.inputs, [model.outputs[0], bias_tensor])

ntest=10
for bno in tqdm(range(ntest)):
    [i0, s0], [o0, bt0] = next(test_gen)
    p0, b0 = bmodel.predict([i0, s0])
    if isinstance(p0, list):
        p0 = p0[0]
    bt.append(bt0[0,...])
    bl.append(b0[0,...])
    pl.append(p0[0,...])
    ol.append(o0[0,...])
    il.append(i0[0,...])
    sl.append(s0[0,...])

inb = np.array(il)
outb = np.array(ol)
pred = np.array(pl)
if pred.shape[-1] > 1:
    pred = pred[...,0:1]

bias = np.array(bl)[...,0]
true_bias = np.array(bt)
seg = np.array(sl)

fv = fs.Freeview(swap_batch_dim=True)
fv.vol(inb, name='synth image', opts=':linked=1:locked=0:auto_adjust_frame_contrast=1')
#fv.vol(seg, name='seg', colormap='lut', opts=':visible=0:locked=1:linked=1:auto_adjust_frame_contrast=1')
fv.vol(outb, name='ground truth', opts=':visible=0:locked=1:linked=1:auto_adjust_frame_contrast=1')
fv.vol(pred, name='corrected', opts=':linked=1:locked=0:visible=1:auto_adjust_frame_contrast=1')
fv.vol(bias, name='bias field', opts=':colormap=heat:heatscale=0,.2:heatscale_offset=1:visible=0:locked=1:linked=1')
fv.vol(true_bias, name='true bias', opts=':colormap=heat:heatscale=0,.2:heatscale_offset=1:visible=0:locked=1:linked=1')
fv.show(title=f'{name}')


if save_model:
    from datetime import date
    fname = f'{hemi}.exvivo.nrounds.{nrounds}.sigma.{sigma}.norm.{str(date.today())}.lsigma.{lsigma}.h5'
    print(f'saving model to {fname}')
    smodel.compile(loss=losses)
    smodel.save(fname)

pfc(write_cb.fname, keys=['loss', 'lr'], close_all=True, smooth=None, 
    remove_outlier_thresh=2, outlier_whalf=4, plot_block=True)

