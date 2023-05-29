#%load_ext autoreload
#%autoreload 2

import os, socket
from tensorflow import keras
import tensorflow as tf
from tensorflow.keras import backend as K
from glob import glob
from tqdm import tqdm

import voxelmorph as vxm
import neurite_sandbox as nes
import neurite as ne
import voxelmorph_sandbox as vxms

use_surfa = True
use_surfa = False
from neurite_sandbox.tf.utils.utils import plot_fit_callback as pfc

if use_surfa:
    import surfa as sf
import fsdeeplearn as fsd
import freesurfer as fs
from fsdeeplearn import surf_utils

# small ones
enc_nf = [32, 32, 32, 32]
dec_nf = [32, 32, 32, 32, 32, 16, 16]
# big ones
enc_nfeats = [64, 96, 128, 128, 128]
dec_nfeats = [128, 128, 96, 64, 64, 32,32]
unet_nfeats = [enc_nfeats, dec_nfeats]
warp_downsize=1


dropout = .0
nfeats = 128
feat_mul = 2
nconv_layers = 30
nfeats = 64
nconv_layers = 15
dilation_rate = 1


nlocal = 1
nfilters = 1
ksize = 1
lc_pad = int((ksize-1)//2)
pad = 16

batch_size = 16

model_dir = 'models'
gpu_id = -1
host = socket.gethostname()
hemi = 'lh'
hemi = 'rh'
atlas_dir='/autofs/space/pac_001/tiamat_copy/3/users/subjects/aparc_atlas'

print(f'host name {socket.gethostname()}')

ngpus = 1 if os.getenv('NGPUS') is None else int(os.getenv('NGPUS'))

print(f'using {ngpus} gpus')
if ngpus > 1:
    model_device = '/gpu:0'
    synth_device = '/gpu:1'
    synth_gpu = 1
    dev_str = "0, 1"
else:
    model_device = '/gpu:0'
    synth_device = '/cpu:0'
    synth_gpu = -1
    dev_str = "0"

os.environ["CUDA_VISIBLE_DEVICES"] = dev_str


fit1=True
fit2=True
fit3=True

fit1=False
fit2=False
fit3=False


print(f'model_device {model_device}, synth_device {synth_device}, dev_str {dev_str}')
print(f'fit1 {fit1}, fit2 {fit2}, fit3 {fit3}, hemi {hemi}')

print(f'physical GPU # is {os.getenv("SLURM_STEP_GPUS")}')
ret = ne.utils.setup_device(dev_str)

#strategy = tf.distribute.MirroredStrategy()
#ngpus = strategy.num_replicas_in_sync
batch_size = (batch_size // ngpus) * ngpus
print('found %d gpus after configuring device, resetting batch_size to %d' % (ngpus, batch_size))

sphere_name='sphere.rot'
curv_name='sulc'
adir = os.path.join(os.getenv('FREESURFER_HOME'), 'average')
fs_atlas_fname = os.path.join(adir, hemi+'.average.curvature.filled.buckner40.tif')
tf.config.get_visible_devices()
if 'inited' not in locals():
    inited = False


if not inited:
    # data loading
    fp = open(os.path.join(atlas_dir, 'aparc-SUBJECTS'))
    subjects = fp.read().split('\n')[:-1]
    fp.close()

    fsdir = os.path.join(os.getenv('FREESURFER_HOME'), 'subjects')
    if use_surfa:   # disabling surfa as it fails
        fs_atlas = sf.load_slice(fs_atlas_fname)
        fsa_surf = sf.load_mesh(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.sphere'))
        fsa_inf = sf.load_mesh(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.inflated'))
    else:
        fs_atlas = fs.Image.read(fs_atlas_fname)
        fsa_surf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.sphere'))
        fsa_inf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.inflated'))


    mrisps_geom = []
    feats = ['inflated.H', 'sulc', 'curv']
    for subject in tqdm(subjects):
        for fno, curv_name in enumerate(feats):
            surf_fname = os.path.join(atlas_dir, subject, 'surf', hemi+'.'+sphere_name)
            curv_fname =os.path.join(atlas_dir, subject, 'surf', hemi+'.'+curv_name)
            mrisp_tmp = surf_utils.loadSphere(surf_fname, curv_fname, padSize=pad)
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
        
    inited = True

lr=1e-4

mrisp_shape = atlas_mean.shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)
warp_sphere_loss = fsd.losses.spherical_loss(tuple(np.array(mrisp_shape_nopad[0:2])//warp_downsize), pad=pad)
sphere_loss = fsd.losses.spherical_loss(mrisp_shape_nopad[0:2], pad=pad, threshold=0, win=[41,41])
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

gen = surf_utils.fsgen(mrisps_geom, atlas_mean, batch_size=batch_size, mean_stream=False)


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

    if use_surfa:
        parcs = [sf.load_overlay(parc_fname) for parc_fname in tqdm(parc_fnames)]
        dkt_spheres = [sf.load_mesh(sphere_fname) for sphere_fname in tqdm(dkt_sphere_fnames)]

        # load the entire dataset to determine labels that exist and compact them
        parcs_all = [sf.load_overlay(parc_fname) for parc_fname in tqdm(parc_fnames_all)]
    else:
        parcs = [fs.Overlay.read(parc_fname) for parc_fname in tqdm(parc_fnames)]
        dkt_spheres = [fs.Surface.read(sphere_fname) for sphere_fname in tqdm(dkt_sphere_fnames)]

        # load the entire dataset to determine labels that exist and compact them
        parcs_all = [fs.Overlay.read(parc_fname) for parc_fname in tqdm(parc_fnames_all)]

    all_labels = nes.py.utils.flatten_lists([list(np.where(parc.data < 0, 0, parc.data)) for parc in parcs_all])
    al = [l if l >= 0 else 0 for l in all_labels]
    lab_to_ind, ind_to_lab = fsd.utils.rebase_labels(al) 

    parcs_training = [lab_to_ind[np.where(parc.data < 0, 0, parc.data)] for parc in parcs]

    # read in and spherically parameterize the geometry and parcellations
    dkt_mrisps_geom = []
    dkt_mrisps_annot = []
    for sno, subject in enumerate(tqdm(dkt_subjects)):
        for fno, curv_name in enumerate(feats): # geometry has 3 scalar fields
            surf_fname = os.path.join(dkt_dir, subject, 'surf', hemi+'.'+sphere_name)
            curv_fname =os.path.join(dkt_dir, subject, 'surf', hemi+'.'+curv_name)
            mrisp_tmp = surf_utils.loadSphere(surf_fname, curv_fname, padSize=pad)
            if fno == 0:
                mrisp = mrisp_tmp[...,np.newaxis]
            else:
                mrisp = np.concatenate((mrisp, mrisp_tmp[...,np.newaxis]), axis=-1)
        dkt_mrisps_geom.append(mrisp)
        if use_surfa:
            mrisp_annot = np.transpose(sf.sphere.SphericalMapNearest(dkt_spheres[sno]).parameterize(parcs_training[sno]),(1,0)).data
        else:
            mrisp_annot = np.transpose(dkt_spheres[sno].parameterize(parcs_training[sno], interp='nearest'),(1,0))
        if pad > 0:
            mrisp_annot = np.pad(mrisp_annot, ((pad,pad),(0,0)),'wrap')
            mrisp_annot = np.pad(mrisp_annot, ((0,0),(pad,pad)),'reflect')

        dkt_mrisps_annot.append(mrisp_annot)

    nclasses = int(np.array(dkt_mrisps_annot).max())+1
    dkt_mrisps_onehot = [fsd.utils.np_one_hot(annot,nclasses) for annot in dkt_mrisps_annot]
    print('using %d classes' % nclasses)
    dkt_mrisps_annot_val = dkt_mrisps_annot[ntraining:ntraining+nval]
    dkt_mrisps_onehot_val = dkt_mrisps_onehot[ntraining:ntraining+nval]
    dkt_mrisps_geom_val = dkt_mrisps_geom[ntraining:ntraining+nval]

    dkt_mrisps_annot = dkt_mrisps_annot[0:ntraining]
    dkt_mrisps_onehot = dkt_mrisps_onehot[0:ntraining]
    dkt_mrisps_geom = dkt_mrisps_geom[0:ntraining]


linear_target = 10
lgen = fsd.surf_utils.fsgen_segreg(dkt_mrisps_geom, atlas_mean, dkt_mrisps_onehot, batch_size=batch_size, mean_stream=False, use_logprob=True)
inbl,outbl = next(lgen)
gen = fsd.surf_utils.fsgen_segreg(dkt_mrisps_geom, atlas_mean, dkt_mrisps_onehot, batch_size=batch_size, mean_stream=False, use_logprob=False)
vbatch = int((1+(nval // ngpus)) * ngpus)
vgen = fsd.surf_utils.fsgen_segreg(dkt_mrisps_geom_val, atlas_mean, dkt_mrisps_onehot_val, batch_size=vbatch, mean_stream=False, use_logprob=False)
lvgen = fsd.surf_utils.fsgen_segreg(dkt_mrisps_geom_val, atlas_mean, dkt_mrisps_onehot_val, batch_size=vbatch, mean_stream=False, use_logprob=True)

print('training segreg atlas with dropout=%2.2f, nfeats=%d and nconvs=%d, dilations=%d' % (dropout, nfeats, nconv_layers, dilation_rate))


warp_smooth_wt = 3
dice_wt = 1
geom_match_wt = 1
which_smooth_loss = 'lap' 
which_smooth_loss = 'grad'

with tf.device(model_device):
    model_semi = vxms.networks.TemplateCreationSemiSupervised(mrisp_shape, nclasses=nclasses, atlas_feats=len(feats), src_feats=len(feats), nb_unet_features=unet_nfeats, mean_cap=100)

    stopping_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=500)
    lr_callback = tf.keras.callbacks.ReduceLROnPlateau(factor=0.8, patience=100, cooldown=10, monitor='loss')
    if which_smooth_loss == 'lap':
        smooth_loss = warp_sphere_loss.laplacianLossEuclidean(penalty='l2')
    else:
        smooth_loss = warp_sphere_loss.gradientLoss(penalty='l2')

    model_semi_name = '%s/%s.fs_atlas.warp%2.2f.geom%2.2f.dice%2.2f.%s.semi.test' % (model_dir, hemi, warp_smooth_wt, geom_match_wt, dice_wt, which_smooth_loss)
    model_semi_fname = model_semi_name + '.h5'
    initial_epoch = 0
    write_cb1 = nes.callbacks.WriteHist(model_semi_name+f'.fit1.txt', mode='w' if initial_epoch == 0 else 'a')
    mc = tf.keras.callbacks.ModelCheckpoint(model_semi_fname, monitor='val_loss', mode='min', save_best_only=True, save_weights_only=True)
    callbacks1 = [stopping_callback, lr_callback, write_cb1]

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

    linear_losses = [channel_loss, l2_loss_inv, l2_loss, smooth_loss]
    linear_loss_weights = [geom_match_wt, geom_match_wt/5, dice_wt/100, warp_smooth_wt] # scale down to account for log-loss

    if fit1:   # do first round of fitting with linear outputs to init things
        print(f'FIT1: saving model weights to {write_cb1.fname}')
        nes.tf.utils.check_and_compile(model_semi.references.model_linear, check_losses=True, optimizer=keras.optimizers.Adam(learning_rate=lr), loss=linear_losses, loss_weights=linear_loss_weights)
        fhist1 = model_semi.references.model_linear.fit(lgen, epochs=30, steps_per_epoch=100, callbacks=callbacks1, validation_data=lvgen, validation_steps=1)
    else:
        if fit2:
            print(f'FIT1: loading model weights from {write_cb1.fname}')
            model_semi.load_weights(mc.filepath)

    model_semi_name = '%s/%s.fs_atlas.warp%2.2f.geom%2.2f.dice%2.2f.%s.semi.test' % (model_dir, hemi, warp_smooth_wt, geom_match_wt, dice_wt, which_smooth_loss)
    model_semi_fname = model_semi_name + '.h5'
    initial_epoch = 0
    write_cb2 = nes.callbacks.WriteHist(model_semi_name+f'.fit2.txt', mode='w' if initial_epoch == 0 else 'a')
    if fit2:  # use softmax as output and dice as loss after init
        losses =        [channel_loss, l2_loss_inv, dice_loss, smooth_loss]
        loss_weights = [geom_match_wt, geom_match_wt/5, dice_wt, warp_smooth_wt]
        callbacks2 = [mc, lr_callback, stopping_callback, write_cb2]
        nes.tf.utils.check_and_compile(model_semi, check_losses=True, optimizer=keras.optimizers.Adam(learning_rate=lr), loss=losses, loss_weights=loss_weights)
        # losses, compile and fit
        print(f'FIT2: saving model weights to {write_cb2.fname}')
        fhist2 = model_semi.fit(gen, epochs=5000, steps_per_epoch=100, callbacks=callbacks2, validation_data=vgen, validation_steps=1)
    else:
        if fit3:
            print(f'FIT2: loading model weights from {mc.filepath}')
            model_semi.load_weights(mc.filepath)

    model_semi_name = '%s/%s.fs_atlas.warp%2.2f.geom%2.2f.dice%2.2f.%s.semi.finetune' % (model_dir, hemi, warp_smooth_wt, geom_match_wt, dice_wt, which_smooth_loss)
    model_semi_fname = model_semi_name + '..h5'
    write_cb3 = nes.callbacks.WriteHist(model_semi_name+f'.fit3..txt', mode='w' if initial_epoch == 0 else 'a')
    if fit3:  # only train atlas
        losses =        [channel_loss, l2_loss_inv, dice_loss, smooth_loss]
        loss_weights = [0, 0, 1, 0]
        mc = tf.keras.callbacks.ModelCheckpoint(model_semi_fname, monitor='val_loss', mode='min', save_best_only=True, save_weights_only=True)
        callbacks3 = [stopping_callback, lr_callback, mc, write_cb3]
        print(f'FIT3: saving model weights to {write_cb3.fname}')

        for layer in model_semi.layers:
            layer.trainable=True if layer.name == 'parc_atlas' else False

        nes.tf.utils.check_and_compile(model_semi, check_losses=True, optimizer=keras.optimizers.Adam(learning_rate=lr), loss=losses, loss_weights=loss_weights)
        fhist3 = model_semi.fit(gen, epochs=5000, steps_per_epoch=100, callbacks=callbacks3, validation_data=vgen, validation_steps=1)
    else:
        print(f'FIT3: loading model weights from {mc.filepath}')
        model_semi.load_weights(mc.filepath)

 
warp_model = tf.keras.Model(model_semi.inputs, model_semi.references.pos_flow)
img_input = tf.keras.Input(shape=dkt_mrisps_onehot[0].shape)
y_img = vxm.layers.SpatialTransformer(interp_method='linear', fill_value=None)([img_input, warp_model.output])
xform_model = tf.keras.Model(warp_model.inputs + [img_input], y_img)
inputs = [np.array(dkt_mrisps_geom), np.array(dkt_mrisps_onehot)]
prior_atlas = xform_model.predict(inputs)
#sf.Slice(np.array(prior_atlas).mean(axis=0)).save(hemi + '.' + 'DKTatlas.priors.mgz')

outputs = model_semi.predict([np.array(dkt_mrisps_geom)])
atlas_mean_in_sub = outputs[0][:,pad:-pad, pad:-pad, :]
sub_in_atlas = outputs[1][:, pad:-pad, pad:-pad, :]
atlas_in_sub_seg = ind_to_lab[np.argmax(outputs[2], axis=-1)[:, pad:-pad, pad:-pad][..., np.newaxis]]
fv = fs.Freeview(swap_batch_dim=True)
fv.vol(np.array(dkt_mrisps_geom)[:, pad:-pad, pad:-pad, :], name='sub geom')
fv.vol(sub_in_atlas, name='sub in atlas geom')
fv.vol(atlas_mean_in_sub, name='atlas geom in sub')
fv.vol(atlas_in_sub_seg, name='atlas seg in sub', opts=':colormap=lut')
fv.show() 

#pfc([write_cb.fname], keys=['val_loss', 'val_atlas_seg_in_sub_lin_loss'],
#    close_all=True, smooth=5, remove_outlier_thresh=2, outlier_whalf=4, plot_block=True)

pfc([write_cb3.fname], keys=['val_loss', 'val_sub_in_atlas_loss', 'val_atlas_in_sub_seg_loss'], 
    close_all=True, smooth=5, remove_outlier_thresh=2, outlier_whalf=4, plot_block=True)
#plt.close('all')
#plt.plot(fh.history['val_loss'])
#plt.plot(fh.history['val_atlas_out_loss'])
#plt.plot(fh.history['val_vxm_dense_transformer_loss'])
#plt.legend(['loss', 'seg loss', 'xform loss'])
#plt.grid()
#plt.show(block=False)
