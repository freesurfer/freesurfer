#%load_ext autoreload
#%autoreload 2

import os, socket
import copy
import scipy

import freesurfer as fs
from freesurfer import deeplearn as fsd
import neurite_sandbox as nes
from freesurfer import deeplearn as fsd

import neurite as ne

from atlasnet import *

batch_size = 16

model_dir = 'models'
gpu_id = -1
host = socket.gethostname()
hemi = 'rh'
atlas_dir='/autofs/space/pac_001/tiamat_copy/3/users/subjects/aparc_atlas'

if os.getenv('NGPUS'):
    ngpus = int(os.getenv('NGPUS'))
    ngpus=1
    gpu_str = '0'
    for g in range(1,ngpus):
        gpu_str += ',%d' % g
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu_str
    print('reading %d GPUS from env and setting CUDA_VISIBLE_DEVICES to %s' % (ngpus, gpu_str))
elif os.getenv('CUDA_VISIBLE_DEVICES'):
    gpu_list = os.getenv('CUDA_VISIBLE_DEVICES')
    ngpus = len(gpu_list.split(','))
elif os.getenv('NGPUS'):
    ngpus = int(os.getenv('NGPUS'))
    gpu_str = '0'
    for g in range(1,ngpus):
        gpu_str += ',%d' % g
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu_str
    print('reading %d GPUS from env and setting CUDA_VISIBLE_DEVICES to %s' % (ngpus, gpu_str))
elif host.startswith('tebo'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3'
elif host.startswith('A100'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3,4,5,6,7'
elif host.startswith('rtx-02'):
    if hemi == 'lh':
        os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2'
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = '3,4,5'
    batch_size=16
elif host.startswith('rtx-0'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3,4,5'
    batch_size=18
elif host.startswith('mlscgpu1'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '2,3,4,5,9'
    batch_size = 15
elif host.startswith('mlscgpu1'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3'
elif host.startswith('A100'):
    os.environ["CUDA_VISIBLE_DEVICES"] = '0,1,2,3'
else:
    gpu_id = 0
    fsd.configure(gpu=gpu_id)

if gpu_id >= 0:
    print('using gpu %d on host %s' % (gpu_id, host))
else:
    gpu_list = os.environ["CUDA_VISIBLE_DEVICES"]
    ngpus = len(gpu_list.split(','))
    batch_size = (batch_size // ngpus) * ngpus
    print('using %d gpus %s on host %s with batch_size %d' % (ngpus,gpu_list, host, batch_size))



from tensorflow import keras
import tensorflow as tf
from tensorflow.keras import backend as K
from surf_utils import *
import socket, os, sys
from glob import glob
from cnn_sphere_register import losspad, prep
import voxelmorph as vxm
import time
import copy

import scipy.io as sio 

strategy = tf.distribute.MirroredStrategy()
ngpus = strategy.num_replicas_in_sync
batch_size = (batch_size // ngpus) * ngpus
print('found %d gpus after configuring device, resetting batch_size to %d' % (ngpus, batch_size))

print('processing %s' % hemi)
sphere_name='sphere.rot'
curv_name='sulc'
adir = os.path.join(os.getenv('FREESURFER_HOME'), 'average')
fs_atlas_fname = os.path.join(adir, hemi+'.average.curvature.filled.buckner40.tif')
tf.config.get_visible_devices()

if 'fs_atlas' not in locals() and 'fs_atlas' not in globals():
    # data loading
    fp = open(os.path.join(atlas_dir, 'aparc-SUBJECTS'))
    subjects = fp.read().split('\n')[:-1]
    fp.close()

    fs_atlas = fs.Image.read(fs_atlas_fname)

    fsdir = os.path.join(os.getenv('FREESURFER_HOME'), 'subjects')
    fsa_surf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.sphere'))
    fsa_inf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.inflated'))


    mrisps_geom = []
    feats = ['inflated.H', 'sulc', 'curv']
    for subject in tqdm(subjects):
        for fno, curv_name in enumerate(feats):
            surf_fname = os.path.join(atlas_dir, subject, 'surf', hemi+'.'+sphere_name)
            curv_fname =os.path.join(atlas_dir, subject, 'surf', hemi+'.'+curv_name)
            mrisp_tmp = prep.spherePad(surf_fname, curv_fname, padSize=pad)
            if fno == 0:
                mrisp = mrisp_tmp[...,np.newaxis]
            else:
                mrisp = np.concatenate((mrisp, mrisp_tmp[...,np.newaxis]), axis=-1)
        mrisps_geom.append(mrisp)


    for fno, curv_name in enumerate(feats):
        atlas_tmp = np.flip(fs_atlas.data[...,fno*3],axis=1).transpose()
        atlas_tmp = (atlas_tmp - np.median(atlas_tmp))/np.std(atlas_tmp)
        atlas_tmp2 = np.sqrt(np.flip(fs_atlas.data[...,fno*3+1],axis=1).transpose())  # convert from variances
        if fno == 0:
            atlas_mean = atlas_tmp[...,np.newaxis]
            atlas_sigma = atlas_tmp2[...,np.newaxis]
        else:
            atlas_mean = np.concatenate((atlas_mean, atlas_tmp[...,np.newaxis]), axis=-1)
            atlas_sigma = np.concatenate((atlas_sigma, atlas_tmp2[...,np.newaxis]), axis=-1)

    if pad > 0:
        atlas_mean = np.pad(atlas_mean, ((pad,pad),(0,0),(0,0)),'wrap')
        atlas_mean = np.pad(atlas_mean, ((0,0),(pad,pad),(0,0)),'reflect')
        atlas_sigma = np.pad(atlas_sigma, ((pad,pad),(0,0),(0,0)),'wrap')
        atlas_sigma = np.pad(atlas_sigma, ((0,0),(pad,pad),(0,0)),'reflect')

lr=1e-4

mrisp_shape = atlas_mean.shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)
warp_sphere_loss = losspad.spherical_loss(tuple(np.array(mrisp_shape_nopad[0:2])//warp_downsize), pad=pad)
sphere_loss = losspad.spherical_loss(mrisp_shape_nopad[0:2], pad=pad, threshold=0, win=[41,41])
dice_loss = sphere_loss.dice_loss(1)

if 0:
    t_pred = tf.zeros((2,)+(mrisp_shape[0:-1]+(2,)))
    lfunc = sphere_loss.laplacianLossEuclidean(penalty='l1')
    l = lfunc(t_pred, t_pred)
    assert 0

if 0:
    t_im = tf.convert_to_tensor(mrisps_geom[0][np.newaxis])
    t_atlas_mean = tf.convert_to_tensor(atlas_mean[np.newaxis])
    t_atlas_sigma = tf.convert_to_tensor(atlas_sigma[np.newaxis])
    l = sphere_loss.l2_loss(1, t_atlas_sigma)(t_atlas_mean, t_atlas_sigma)

gen = fsgen(mrisps_geom, atlas_mean, batch_size=batch_size, mean_stream=False)



#with strategy.scope():


dkt_dir = os.path.join(atlas_dir, 'DKTatlas')
parc_name = 'labels.DKT31.manual.annot'
parc_name = 'labels.DKT31.manual.2.annot'
parc_fnames_all = glob(os.path.join(dkt_dir, '*/label', hemi+'.'+parc_name))
ntraining = 50
nval = 10
parc_fnames = parc_fnames_all[:ntraining+nval]

if 'read_hemi' not in locals() and 'read_hemi' not in globals():
    read_hemi = hemi

# check to see if data was loaded, and if not load it in
if read_hemi is not hemi or ('dkt_mrisps_onehot' not in locals() and 'dkt_mrisps_onehot' not in globals()):
    # load in all the surfaces and parcellations for training data
    dkt_subjects = [parc_fname.split('/')[-3] for parc_fname in parc_fnames]
    dkt_sphere_fnames = [os.path.join(dkt_dir, subject, 'surf', hemi+'.'+sphere_name) for subject in dkt_subjects]

    parcs = [fs.Overlay.read(parc_fname) for parc_fname in tqdm(parc_fnames)]
    dkt_spheres = [fs.Surface.read(sphere_fname) for sphere_fname in tqdm(dkt_sphere_fnames)]

    # load the entire dataset to determine labels that exist and compact them
    parcs_all = [fs.Overlay.read(parc_fname) for parc_fname in tqdm(parc_fnames_all)]
    all_labels = nes.py.utils.flatten_lists([list(np.where(parc.data<0, 0, parc.data)) for parc in parcs_all])
    al = [l if l >= 0 else 0 for l in all_labels]
    lab_to_ind, ind_to_lab = rebase_labels(al) 

    parcs_training = [lab_to_ind[np.where(parc.data<0, 0, parc.data)] for parc in parcs]

    # read in and spherically parameterize the geometry and parcellations
    dkt_mrisps_geom = []
    dkt_mrisps_annot = []
    for sno, subject in enumerate(tqdm(dkt_subjects)):
        for fno, curv_name in enumerate(feats): # geometry has 3 scalar fields
            surf_fname = os.path.join(dkt_dir, subject, 'surf', hemi+'.'+sphere_name)
            curv_fname =os.path.join(dkt_dir, subject, 'surf', hemi+'.'+curv_name)
            mrisp_tmp = prep.spherePad(surf_fname, curv_fname, padSize=pad)
            if fno == 0:
                mrisp = mrisp_tmp[...,np.newaxis]
            else:
                mrisp = np.concatenate((mrisp, mrisp_tmp[...,np.newaxis]), axis=-1)
        dkt_mrisps_geom.append(mrisp)
        mrisp_annot = np.transpose(dkt_spheres[sno].parameterize(parcs_training[sno], interp='nearest'),(1,0))
        if pad > 0:
            mrisp_annot = np.pad(mrisp_annot, ((pad,pad),(0,0)),'wrap')
            mrisp_annot = np.pad(mrisp_annot, ((0,0),(pad,pad)),'reflect')

        dkt_mrisps_annot.append(mrisp_annot)

    nclasses = int(np.array(dkt_mrisps_annot).max())+1
    dkt_mrisps_onehot = [np_one_hot(annot,nclasses) for annot in dkt_mrisps_annot]
    print('using %d classes' % nclasses)
    dkt_mrisps_annot_val = dkt_mrisps_annot[ntraining:ntraining+nval]
    dkt_mrisps_onehot_val = dkt_mrisps_onehot[ntraining:ntraining+nval]
    dkt_mrisps_geom_val = dkt_mrisps_geom[ntraining:ntraining+nval]

    dkt_mrisps_annot = dkt_mrisps_annot[0:ntraining]
    dkt_mrisps_onehot = dkt_mrisps_onehot[0:ntraining]
    dkt_mrisps_geom = dkt_mrisps_geom[0:ntraining]


linear_target = 10
lgen = fsgen_segreg(dkt_mrisps_geom, atlas_mean, dkt_mrisps_onehot, batch_size=batch_size, mean_stream=False, use_logprob=True)
inbl,outbl = next(lgen)
gen = fsgen_segreg(dkt_mrisps_geom, atlas_mean, dkt_mrisps_onehot, batch_size=batch_size, mean_stream=False, use_logprob=False)
vbatch = int((1+(nval // ngpus)) * ngpus)
vgen = fsgen_segreg(dkt_mrisps_geom_val, atlas_mean, dkt_mrisps_onehot_val, batch_size=vbatch, mean_stream=False, use_logprob=False)
lvgen = fsgen_segreg(dkt_mrisps_geom_val, atlas_mean, dkt_mrisps_onehot_val, batch_size=vbatch, mean_stream=False, use_logprob=True)

print('training segreg atlas with dropout=%2.2f, nfeats=%d and nconvs=%d, dilations=%d' % (dropout, nfeats, nconv_layers, dilation_rate))


warp_smooth_wt = 3
dice_wt = 1
geom_match_wt = 1
which_smooth_loss = 'lap' 
which_smooth_loss = 'grad'

fit1=True
fit2=True
fit3=True

import tensorflow.keras.initializers as KI
with strategy.scope():
    model_semi = vxms.networks.TemplateCreationSemiSupervised(mrisp_shape, nclasses=nclasses, atlas_feats=len(feats), src_feats=len(feats), nb_unet_features=unet_nfeats, mean_cap=100)

    stopping_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=500)
    lr_callback = tf.keras.callbacks.ReduceLROnPlateau(factor=0.8, patience=100, cooldown=10, monitor='loss')
    if which_smooth_loss == 'lap':
        smooth_loss = warp_sphere_loss.laplacianLossEuclidean(penalty='l2')
    else:
        smooth_loss = warp_sphere_loss.gradientLoss(penalty='l2')

    model_semi_fname = '%s/%s.fs_atlas.warp%2.2f.geom%2.2f.dice%2.2f.%s.semi.test.h5' % (model_dir, hemi, warp_smooth_wt, geom_match_wt, dice_wt, which_smooth_loss)
    mc = tf.keras.callbacks.ModelCheckpoint(model_semi_fname, monitor='val_loss', mode='min', save_best_only=True, save_weights_only=True)
    callbacks = [stopping_callback, lr_callback]
    print('saving model weights to %s' % model_semi_fname)

    # create a likelihood loss for all the geometric measures
    # and split it into different channels so that they can be weighted
    # differently. Inverse losses are just mse with no variance estimates
    l2_loss_fwd = sphere_loss.l2_loss(1,atlas_sigma)
    l2_loss_c0 = sphere_loss.l2_loss(1,atlas_sigma[...,0:1])
    l2_loss_c1 = sphere_loss.l2_loss(1,atlas_sigma[...,1:2])
    l2_loss_c2 = sphere_loss.l2_loss(1,atlas_sigma[...,2:])

    l2_loss_inv = sphere_loss.l2_loss(1,None)
    l2_loss = sphere_loss.l2_loss(nclasses,None)

    # weight sulc in frame 1 more heavily
    channel_loss = nes.losses.channelwise_losses_decorator([0,1,2], [l2_loss_c0, l2_loss_c1, l2_loss_c2], loss_weights=[.1,1,.1],norm_weights=True)

    linear_losses = [channel_loss, l2_loss_inv, smooth_loss, l2_loss]
    linear_loss_weights = [geom_match_wt, geom_match_wt/5, warp_smooth_wt, 1.0/100] # scale down to account for log-loss

    if fit1:   # do first round of fitting with linear outputs to init things
        nes.tf.utils.check_and_compile(model_semi.references.model_linear, check_losses=True, optimizer=keras.optimizers.Adam(lr=lr), loss=linear_losses, loss_weights=linear_loss_weights)
        fhist1 = model_semi.references.model_linear.fit(lgen, epochs=30, steps_per_epoch=100, callbacks=callbacks, validation_data=lvgen, validation_steps=1)

    if fit2:  # use softmax as output and dice as loss after init
        losses =        [channel_loss, l2_loss_inv, smooth_loss, dice_loss]
        loss_weights = [geom_match_wt, geom_match_wt/5, warp_smooth_wt, dice_wt]
        callbacks = [mc, lr_callback, stopping_callback]
        nes.tf.utils.check_and_compile(model_semi, check_losses=True, optimizer=keras.optimizers.Adam(lr=lr), loss=losses, loss_weights=loss_weights)
        # losses, compile and fit
        fhist2 = model_semi.fit(gen, epochs=5000, steps_per_epoch=100, callbacks=callbacks, validation_data=vgen, validation_steps=1)
    else:
        fhist2 = []
        model_semi.load_weights(model_semi_fname)

    if fit3:  # only train atlas
        losses =        [channel_loss, l2_loss_inv, smooth_loss, dice_loss]
        loss_weights = [0, 0, 0, 1]
        model_semi_fname = '%s/%s.fs_atlas.warp%2.2f.geom%2.2f.dice%2.2f.%s.semi.finetune.h5' % (model_dir, hemi, warp_smooth_wt, geom_match_wt, dice_wt, which_smooth_loss)
        mc = tf.keras.callbacks.ModelCheckpoint(model_semi_fname, monitor='val_loss', mode='min', save_best_only=True, save_weights_only=True)
        callbacks = [stopping_callback, lr_callback, mc]
        print('saving model weights to %s' % model_semi_fname)

        for layer in model_semi.layers:
            layer.trainable=True if layer.name == 'parc_atlas' else False

        nes.tf.utils.check_and_compile(model_semi, check_losses=True, optimizer=keras.optimizers.Adam(lr=lr), loss=losses, loss_weights=loss_weights)
        fhist3 = model_semi.fit(gen, epochs=5000, steps_per_epoch=100, callbacks=callbacks, validation_data=vgen, validation_steps=1)


 
warp_model = tf.keras.Model(model_semi.inputs, model_semi.references.pos_flow)
img_input = tf.keras.Input(shape=dkt_mrisps_onehot[0].shape)
y_img = vxm.layers.SpatialTransformer(interp_method='linear', fill_value=None)([img_input, warp_model.output])
xform_model = tf.keras.Model(warp_model.inputs + [img_input], y_img)
inputs = [np.array(dkt_mrisps_geom), np.array(dkt_mrisps_onehot)]
prior_atlas = xform_model.predict(inputs)
fs.Image(np.array(prior_atlas).mean(axis=0)).write(hemi+'.'+'DKTatlas.priors.mgz')

if fit3:
    fh = fhist3
else:
    fh = fhist2

plt.close('all')
plt.plot(fh.history['val_loss'])
plt.plot(fh.history['val_atlas_out_loss'])
plt.plot(fh.history['val_vxm_dense_transformer_loss'])
plt.legend(['loss', 'seg loss', 'xform loss'])
plt.grid()
plt.show(block=False)


