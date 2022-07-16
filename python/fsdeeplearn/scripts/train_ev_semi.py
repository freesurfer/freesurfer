#%load_ext autoreload
#%autoreload 2

import os, socket, sys, copy
from glob import glob
from tensorflow import keras
from tensorflow.keras import backend as K
from tensorflow.keras import layers as KL
import tensorflow as tf

from tqdm import tqdm

import freesurfer as fs
from freesurfer import deeplearn as fsd
import tensorflow as tf
import neurite_sandbox as nes
import neurite as ne
import voxelmorph as vxm

from neurite_sandbox.tf.utils.utils import plot_fit_callback as pfc

npos_encodings=16
npos_encodings=0
npos_encodings=4

save_images = True
save_images = False

dofit = True
dofit = False

which_training = 'func'
which_training = 'geom'

which_fmri = 'POW'
which_fmri = 'L2'
which_fmri = 'NCC'

which_geom = 'POW'
which_geom = 'NCC'
which_geom = 'L2'

which_atlas = 'NCC'

batch_size = 8

warp_smooth_wt = 1
warp_smooth_wt = .1
geom_match_wt = .1
geom_match_wt = 1

which_target = 'subject'
which_target = 'atlas'

hemi = 'lh'

build_atlases = True
build_atlases = False

smooth_steps = 200  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))
smooth_steps = 50  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))
smooth_steps = 1  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))

func_wt_list = [.1]
func_wt_list = [100]
func_wt_list = [10]
func_wt_list = [1]
func_wt_list = [0]
func_wt_list = [5]
func_wt_list = [90]

model_dir = 'models'
name = f'{model_dir}/{hemi}.surf2surf.{which_training}.warp_smooth_wt_%2.2f.geom_wt_{geom_match_wt}.smooth_{smooth_steps}.geom_{which_geom}.fmri_{which_fmri}' % (warp_smooth_wt)
aname = f'{model_dir}/{hemi}.surf2surf.semi.warp_smooth_wt_%2.2f.geom_wt_{geom_match_wt}.smooth_{smooth_steps}.{which_fmri}.npos_{npos_encodings}' % (warp_smooth_wt)

gpuid = -1
host = socket.gethostname()
if os.getenv('NGPUS'):
    ngpus = int(os.getenv('NGPUS'))
#    ngpus=1
    gpu_str = '0'
    for g in range(1,ngpus):
        gpu_str += ',%d' % g
    if ngpus == 1:
        gpuid = 0
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu_str
    print('reading %d GPUS from env and setting CUDA_VISIBLE_DEVICES to %s' % (ngpus, gpu_str))
elif os.getenv('CUDA_VISIBLE_DEVICES'):
    gpu_list = os.getenv('CUDA_VISIBLE_DEVICES')
    ngpus = len(gpu_list.split(','))
else:
    gpuid = 0
    fsd.configure(gpu=gpuid)

if gpuid >= 0:
    print('using gpu %d on host %s' % (gpuid, host))
else:
    gpu_list = os.environ["CUDA_VISIBLE_DEVICES"]
    ngpus = len(gpu_list.split(','))
    batch_size = (batch_size // ngpus) * ngpus
    print('using %d gpus %s on host %s with batch_size %d' % (ngpus,gpu_list, host, batch_size))

print(f'dofit {dofit}, which_geom {which_geom}, which_fmri {which_fmri}, which_training {which_training}, smooth_steps {smooth_steps}, target {which_target}, build_atlases {build_atlases}')

devices = tf.config.get_visible_devices()
# device, ngpus = ne.utils.setup_device(gpuid=gpuid)
# tf.compat.v1.disable_eager_execution()

# paths to data and list of subjects and such
base_dir = '/autofs/cluster/p41/resting_state/Freesurfer_Surface'
sdir1 = 'CG1'
sdir2 = 'CG2'
fsdir = os.path.join(base_dir, sdir1, 'FS')
fsdir_val = os.path.join(base_dir, sdir2, 'FS')

sphere_name = '%s.sphere.rot' % hemi
sulc_name = '%s.sulc' % hemi

fsa_sphere = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.sphere.reg'))
fsa_inf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.inflated'))

contrast_name = 'SvsN'
func_dir = os.path.join('bold','bold.self.sm0.%s.lang' % hemi, contrast_name, 't.nii.gz')
spaths = glob(os.path.join(fsdir, 'sub???'))
crop = -1
crop2 = -1
spaths = spaths[0:crop]
spaths_val = glob(os.path.join(fsdir_val, 'sub???'))[0:crop2]

 # load all the data (if not previously loaded)
pad = 8
smooth_read = -1

geom_names = ['inflated.H', 'sulc', 'curv']
def get_angles(pos, i, d_model):
  angle_rates = 1 / np.power(10000, (2 * (i//2)) / np.float32(d_model))
  return pos * angle_rates

def pos_encoding2D(inshape, npos):
    pe = np.zeros(inshape + (2*npos,))

    for pno in range(npos):
        wavelen = inshape[0] / (pno+1)
        const = 2 * np.pi / wavelen
        pex = np.repeat(np.sin(np.arange(inshape[0]) * const)[...,np.newaxis], inshape[1], axis=1)
        pe[..., pno] = pex
        wavelen = inshape[1] / (pno+1)
        const = 2 * np.pi / wavelen
        pey = np.repeat(np.sin(np.arange(inshape[1]) * const)[np.newaxis], inshape[0], axis=0)
        pe[..., npos+pno] = pey
                                
    return pe


if 'spheres' not in locals() and 'spheres' not in globals() and smooth_steps != smooth_read:
    geom_fname = f'mrisps_geom_smooth_{smooth_steps}.npy'
    func_fname = f'mrisps_func_smooth_{smooth_steps}.npy'
    # func_fname = f'mrisps_func_smooth_50.npy'
    spheres_fname = f'spheres_smooth_{smooth_steps}.npy'
    snames_fname = f'snames_smooth_{smooth_steps}.npy'
    geom_fname_val = f'mrisps_geom_smooth_{smooth_steps}_val.npy'
    func_fname_val = f'mrisps_func_smooth_{smooth_steps}_val.npy'
    # func_fname_val = f'mrisps_func_smooth_50.npy'
    spheres_fname_val = f'spheres_smooth_{smooth_steps}_val.npy'
    snames_fname_val = f'snames_smooth_{smooth_steps}_val.npy'
    if os.path.exists(geom_fname):
        print('reading lists from npy cache...')
        mrisps_geom = np.load(geom_fname)
        mrisps_func = np.load(func_fname)
        spheres = np.load(spheres_fname, allow_pickle=True)
        snames = np.load(snames_fname)

        mrisps_geom_val = np.load(geom_fname_val)
        mrisps_func_val = np.load(func_fname_val)
        spheres_val = np.load(spheres_fname_val, allow_pickle=True)
        snames_val = np.load(snames_fname_val)
    else:
        mrisps_geom, mrisps_func, spheres, snames = fsd.surf_utils.load_func_and_spheres(spaths, base_dir, func_dir, sdir1, 
                                                                                         sphere_name, hemi, pad=pad, 
                                                                                         smooth_steps=smooth_steps,
                                                                                         geom_names=geom_names)
        mrisps_geom_val, mrisps_func_val, spheres_val, snames_val = fsd.surf_utils.load_func_and_spheres(spaths_val, base_dir, 
                                                                                                         func_dir, sdir2, 
                                                                                                         sphere_name, hemi, 
                                                                                                         pad=pad, 
                                                                                                         smooth_steps=smooth_steps,
                                                                                                         geom_names=geom_names)
        print('saving lists to npy cache...')
        np.save(geom_fname, mrisps_geom, allow_pickle=True)
        np.save(func_fname, mrisps_func, allow_pickle=True)
        np.save(spheres_fname, spheres, allow_pickle=True)
        np.save(snames_fname, snames, allow_pickle=True)

        np.save(geom_fname_val, mrisps_geom_val, allow_pickle=True)
        np.save(func_fname_val, mrisps_func_val, allow_pickle=True)
        np.save(spheres_fname_val, spheres_val, allow_pickle=True)
        np.save(snames_fname_val, snames_val, allow_pickle=True)
                              
    smooth_read = smooth_steps

nsubjects = len(spaths)

mrisps_func_orig = mrisps_func

# optional - threshold the functional data to remove noise
threshold = 0
if threshold > 0:
    mrisps_func = []
    for mrisp in mrisps_func_orig:
        ind_off = np.nonzero(np.abs(mrisp) < threshold)
        mrisp = copy.deepcopy(mrisp)
        mrisp[ind_off] = 0
        mrisps_func.append(mrisp)


warp_downsize=2
if which_target == 'atlas':
    fname_folding = hemi + f'.{which_atlas}.atlas.folding'
    fname_fmri = hemi + f'.{which_atlas}.atlas.fmri'
    atlas_folding_geom = fs.Image.read(fname_folding + '.mean_geom.mrisp.mgz')
    atlas_folding_func = fs.Image.read(fname_folding + '.mean_func.mrisp.mgz')
    atlas_folding_geom_std = fs.Image.read(fname_folding + '.std_geom.mrisp.mgz')
    atlas_folding_func_std = fs.Image.read(fname_folding + '.std_func.mrisp.mgz')

    atlas_fmri_geom = fs.Image.read(fname_fmri + '.mean_geom.mrisp.mgz')
    atlas_fmri_func = fs.Image.read(fname_fmri + '.mean_func.mrisp.mgz')
    atlas_fmri_geom_std = fs.Image.read(fname_fmri + '.std_geom.mrisp.mgz')
    atlas_fmri_func_std = fs.Image.read(fname_fmri + '.std_func.mrisp.mgz')
    atlas_fmri_func_std.data[atlas_fmri_func_std.data < .01] = .01
    if which_fmri == 'NCC' and 0:
        atlas_fmri_func_norm = atlas_fmri_func.data / atlas_fmri_func_std.data
    else:
        atlas_fmri_func_norm = atlas_fmri_func.data

    fmri_thresh = 2
    fmri_thresh = 0
    fmri_spatial_weights = atlas_fmri_func_norm.copy()
    fmri_spatial_weights[atlas_fmri_func_norm < fmri_thresh] = 0
    fmri_spatial_weights = None

    noise_max = .8
    fmri_noise = 2
    fmri_noise = 4
    fmri_noise = 8
    aug_types = ['noise', 'translate']
    aug_types = []
    aug_types = ['noise']
    # use atlases built from geometry-only
    fgen_fold = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_func, mrisps_geom, atlas_folding_func.data[...,np.newaxis], 
                                                    atlas_folding_geom.data, batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize)
    ggen_fold = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, atlas_folding_geom.data, 
                                                    atlas_folding_func.data[...,np.newaxis], batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize, aug_types=aug_types,
                                                    noise_max=noise_max)
    vgen_fold = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom_val, mrisps_func_val, atlas_folding_geom.data, 
                                                    atlas_folding_func.data[...,np.newaxis], batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize)
    # use atlases built from fmri data also
    fgen_fmri = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_func, mrisps_geom, atlas_fmri_func_norm[...,np.newaxis], 
                                                    atlas_fmri_geom.data, batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize)
    ggen_fmri = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, atlas_fmri_geom.data, 
                                                    atlas_fmri_func_norm[...,np.newaxis], batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize, aug_types=aug_types,
                                                    noise_max=noise_max, fmri_noise=fmri_noise)
    vgen_fmri = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom_val, mrisps_func_val, atlas_fmri_geom.data, 
                                                    atlas_fmri_func_norm[...,np.newaxis], batch_size=batch_size, 
                                                    use_rand=True, warp_downsize=warp_downsize)
else:
    atlas = None

    fgen = fsd.surf_utils.mrisp_semi_gen(mrisps_func, mrisps_geom, batch_size=batch_size, use_rand=True,
                                         warp_downsize=warp_downsize)
    ggen = fsd.surf_utils.mrisp_semi_gen(mrisps_geom, mrisps_func, batch_size=batch_size, use_rand=True, 
                                         warp_downsize=warp_downsize)


if which_training == 'geom':
    ninputs = mrisps_geom[0].shape[-1] 
    semi_inputs = 1
else:
    semi_inputs = mrisps_geom[0].shape[-1] 
    ninputs = 1
    
mrisp_shape = mrisps_geom[0].shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)
warp_sphere_loss = fsd.losses.spherical_loss(tuple(np.array(mrisp_shape_nopad)//warp_downsize), pad=pad//2)
sphere_loss = fsd.losses.spherical_loss(mrisp_shape_nopad, pad=pad, threshold=threshold, win=[31,31])
if which_fmri == 'NCC':
    fmri_loss = sphere_loss.NCC_signed_loss(weight=1, weights=fmri_spatial_weights)
elif which_fmri == 'L2':
    atlas_func_std = atlas_fmri_func_std.data if func_wt_list[0] > 0 else atlas_folding_func_std.data
    fmri_loss = sphere_loss.l2_loss(1, atlas_func_std[...,np.newaxis])
else:
    fmri_loss = sphere_loss.power_loss

if which_geom == 'NCC':
    geom_loss = sphere_loss.NCC_signed_loss(weight=1)
elif which_geom == 'L2':
    atlas_std = atlas_fmri_geom_std.data if func_wt_list[0] > 0 else atlas_folding_geom_std.data
    inflated_H_wt = .01
    curv_wt = 1
    atlas_std[..., 0] *=  (1.0 / inflated_H_wt)
    atlas_std[..., 2] *=  (1.0 / curv_wt)
    geom_loss = sphere_loss.l2_loss(1, atlas_std)
else:
    geom_loss = sphere_loss.power_loss

fhists = []

# vxm unet architecture
fscale = 4
enc_nfeats = [64, 96, 128, 128]
dec_nfeats = [128, 64, 64, 64, 32,32]
unet_nfeats = [[nf * fscale for nf in enc_nfeats], [nf * fscale for nf in dec_nfeats]]

lr = 1e-4
models = []
if not dofit:
    model_names = []
    fit_names = []
    if which_training == 'func':
        func_wt_list = [100]
    else:
        func_wt_list = [5]
        func_wt_list = [0,.1,1,5, 10,100]
        func_wt_list = [0, 1, 5, 10,90, 100]

for func_wt in func_wt_list:
    if func_wt > 1:
        func_wt_this_round = 1
        geom_match_wt_this_round = geom_match_wt / func_wt
    else:
        func_wt_this_round = func_wt
        geom_match_wt_this_round = 1

    strategy = tf.distribute.MirroredStrategy()

    name_this_round = name + f'.func_wt_%2.2f' % func_wt

    # callbacks
    write_cb = nes.callbacks.WriteHist(name_this_round+f'.txt', mode='w')
    lr_cb = nes.tf.callbacks.ReduceLRWithModelCheckpointAndRecovery(fname=name_this_round+'.h5', save_weights_only=False, 
        factor=0.8, patience=100, cooldown=10, monitor='loss', burn_in=50, thresh_inc_factor=3, recovery_decrease_factor=1)
    # with strategy.scope():
    if 1:
        source_input = KL.Input(mrisp_shape+(ninputs,), name='geom_input')
        atlas_input = KL.Input(mrisp_shape+(ninputs,), name='atlas_input')
        fmri_input = KL.Input(mrisp_shape+(semi_inputs,), name='fmri_input')
        source_normed = KL.BatchNormalization(name='source_norm')(source_input)
        fmri_normed = KL.BatchNormalization(name='fmri_norm')(fmri_input)
        input_model = tf.keras.Model([source_input, atlas_input], [source_normed, atlas_input])
        input_model = None
        if npos_encodings > 0:
            model_vxm = vxm.networks.VxmDenseSemiSupervisedSeg(
                mrisp_shape, 
                semi_inputs, 
                src_feats=ninputs,
                trg_feats=ninputs,
                nb_unet_features=unet_nfeats, 
                input_model=input_model,
                int_resolution=warp_downsize, 
                pos_encoding=npos_encodings,
                seg_resolution=1)
        else:
            model_vxm = vxm.networks.VxmDenseSemiSupervisedSeg(
                mrisp_shape, 
                semi_inputs, 
                src_feats=ninputs,
                trg_feats=ninputs,
                input_model=input_model,
                nb_unet_features=unet_nfeats, 
                int_resolution=warp_downsize, 
                seg_resolution=1)

        if dofit:
            print('training with %s loss weight %2.2f and which_fmri = %s, geom_wt=%2.2f, warp_smooth_wt=%2.2f' % 
                  (which_training, func_wt_this_round,which_fmri, geom_match_wt_this_round,warp_smooth_wt))

            new_out1 = KL.Lambda(lambda x : x, name='warp')(model_vxm.outputs[1])
            if which_training == 'func':
                new_out0 = KL.Lambda(lambda x : x, name='fmri')(model_vxm.outputs[0])
                new_out2 = KL.Lambda(lambda x : x, name='geom')(model_vxm.outputs[2])
                losses = [ fmri_loss, warp_sphere_loss.gradientLoss(penalty='l2'), geom_loss]
            else:
                new_out0 = KL.Lambda(lambda x : x, name='geom')(model_vxm.outputs[0])
                new_out2 = KL.Lambda(lambda x : x, name='fmri')(model_vxm.outputs[2])
                losses = [ geom_loss, warp_sphere_loss.gradientLoss(penalty='l2'), fmri_loss]

            model_semi = tf.keras.Model(model_vxm.inputs, [new_out0, new_out1, new_out2])
            model_semi.references = model_vxm.references
            # losses and compile
            if which_training == 'func':
                loss_weights = [func_wt_this_round,warp_smooth_wt,geom_match_wt_this_round]
            else:
                loss_weights = [geom_match_wt_this_round,warp_smooth_wt,func_wt_this_round]

            if which_training == 'func': # function goes through convs
                gen = fgen_fmri if func_wt_this_round > .001 or which_target == 'subject' else fgen_fold
                vgen = fgen_fmri
            else:
                gen = ggen_fmri if func_wt_this_round > .001 or which_target == 'subject' else ggen_fold
                vgen = vgen_fmri

            nes.utils.check_and_compile(model_semi, gen, optimizer=keras.optimizers.Adam(lr=lr), 
                                        loss=losses, loss_weights=loss_weights)
            print(f'writing fit history to {write_cb.fname} and model weights to {lr_cb.fname}')
            callbacks = [write_cb, lr_cb]
            fhist = model_semi.fit(gen, epochs=10000, steps_per_epoch=200, callbacks=callbacks, 
                                   validation_data=vgen, validation_steps=8)
            # fhist = model_semi.fit(gen, epochs=10000, steps_per_epoch=200, callbacks=callbacks, initial_epoch=5000)

            fhists.append(fhist)
            # model_semi.save_weights(model_fname)
        else:
            model_semi = model_vxm
            model_names.append(lr_cb.fname)
            models.append(model_semi)
            fit_names.append(write_cb.fname)
            #pfc(write_cb.fname, keys=['loss', 'geom_loss', 'fmri_loss'], close_all=True, 
            #    smooth=None,  remove_outlier_thresh=2, outlier_whalf=4, plot_block=False)

            if 0:
                model_semi.load_weights(lr_cb.fname)
                inb, outb = next(gen)
                pred = model_semi.predict(inb)
                moving = np.transpose(inb[0], (1,2,3,0))
                fixed = np.transpose(inb[1], (1,2,3,0))
                moving_func = np.transpose(inb[2], (1,2,3,0))
                fixed_func = np.transpose(outb[2], (1,2,3,0))
                moving_warped = np.transpose(pred[0], (1,2,3,0))
                moving_func_warped = np.transpose(pred[2], (1,2,3,0))
                fv = fs.Freeview()
                fv.vol(fixed, name='fixed', opts=':locked=1')
                fv.vol(moving, name='moving', opts=':locked=1:visible=0')
                fv.vol(moving_warped, name='moving_warped', opts=':locked=1')
                fv.vol(moving_func, name='func_moving', colormap='heat', opts=':locked=1:visible=0:linked=1')
                fv.vol(moving_func_warped, name='moving_func_warped_%2.2f' % func_wt, colormap='heat', opts=':locked=1:visible=0:linked=1')
                fv.vol(fixed_func, name='func_fixed', colormap='heat', opts=':locked=1:visible=0:linked=1')
                fv.vol((moving_func+fixed_func)/2, name='avg_b4', colormap='heat', opts=':locked=0:visible=0:linked=1')
                fv.vol((moving_func_warped+fixed_func)/2, name='avg_after', colormap='heat', opts=':locked=0:visible=1:linked=1')
                fv.show(title='func_wt=%2.1f' % func_wt)


if build_atlases:
    fname_geom = 'models/lh.surf2surf.semi.warp_smooth_wt_0.10.geom_wt_1.smooth_1.NCC.func_wt_0.00.h5'
    model_semi.load_weights(fname_geom)

    fname_func = 'models/lh.surf2surf.sup.warp_smooth_wt_0.10.geom_wt_1.smooth_1.NCC.func_wt_10.00.h5'
    model_func = vxm.networks.VxmDenseSemiSupervisedSeg(
        mrisp_shape, 
        mrisps_geom[0].shape[-1], 
        src_feats=1,
        trg_feats=1,
        nb_unet_features=unet_nfeats, 
        int_resolution=warp_downsize, 
        seg_resolution=1)
        
    model_func.load_weights(fname_func)
    func_wt_list = [0, 10]
    model_geom = model_semi
    model_names = [fname_geom, fname_func]
    fit_names = [".".join(fname_geom.split('.')[0:-1])+'.txt', ".".join(fname_func.split('.')[0:-1])+'.txt']


doval = False
doval = True
if doval:
    mrisps_func_test = mrisps_func_val
    snames_test = snames_val
    mrisps_geom_test = mrisps_geom_val
else:
    mrisps_func_test = mrisps_func
    snames_test = snames
    mrisps_geom_test = mrisps_geom

if not dofit:
    mean_geoms = []
    std_geoms = []
    mean_funcs = []
    std_funcs = []
    geom_list = []
    func_list = []
    for fno, func_wt in enumerate(func_wt_list):    # for every different model
        if build_atlases:
            model_semi = model_func if func_wt > 0 else model_geom

        model_semi = models[fno]
        model_semi.load_weights(model_names[fno])
        warped_funcs = []
        warped_geoms = []
        avg_warps = []
        avg_geoms = []
        avg_funcs = []
        trans_geom_model = vxm.networks.Transform(mrisp_shape, nb_feats=mrisps_geom[0].shape[-1])
        trans_func_model = vxm.networks.Transform(mrisp_shape)
        registration_model = model_semi.get_registration_model()
        if which_target == 'subject':
            print(f'computing subject warp matrix for func_wt {func_wt}: {fno+1} of {len(func_wt_list)}')
        else:
            print(f'warping each subject to the atlas for func_wt {func_wt}: {fno+1} of {len(func_wt_list)}')
        for sno, sname in enumerate(tqdm(snames_test)):   # compute and apply warp to implicit midspace for each subject
            if which_target == 'subject':
                wgeoms = []
                inputs = [[],[],[]]
                for sno2, sname2 in enumerate(snames_test):   # compute warp from subject sno to every other subject, then avg it
                    if which_training == 'func':
                        inputs[0].append(mrisps_func_test[sno])
                        inputs[1].append(mrisps_func_test[sno2])
                        inputs[2].append(mrisps_geom_test[sno])
                    else:
                        inputs[0].append(mrisps_geom_test[sno])
                        inputs[1].append(mrisps_geom_test[sno2])
                        inputs[2].append(mrisps_func_test[sno])

                # with strategy.scope():
                if 1:
                    model_inputs = [np.array(inputs[0]), np.array(inputs[1]), np.array(inputs[2])]
                    pred = model_semi.predict(model_inputs)
                
                    if which_training == 'func':
                        wgeoms = pred[2]
                        wfuncs = pred[0]
                    else:
                        wgeoms = pred[0]
                        wfuncs = pred[2]

                    warps = pred[1]
                    warps = registration_model.predict(model_inputs[0:2])
                
                    warped_funcs.append(wfuncs)  # list of all warped funcs for this subject
                    warped_geoms.append(wgeoms)  # list of all warped geoms for this subject

                    # compute the warp to the implicit common space for this subject
                    avg_warp = np.array(warps).mean(axis=0).squeeze()  # this warp takes this subject to the avg space
                    avg_geom = np.array(wgeoms).mean(axis=0).squeeze() 
                    avg_func = np.array(wfuncs).mean(axis=0).squeeze()
                    avg_warps.append(avg_warp)

                    # warp this subject;s geometry and functional data to common space
                    warped_geom = trans_geom_model.predict([mrisps_geom[sno][np.newaxis], avg_warp[np.newaxis]])
                    warped_func = trans_func_model.predict([mrisps_func[sno][np.newaxis], avg_warp[np.newaxis]])
                    
                    avg_geoms.append(warped_geom)    # this subject's geom mapped to average space
                    avg_funcs.append(warped_func)    # this subject's func mapped to average space
        
                # end subject loop
            else:    # target is atlas. Only need to warp each subject once
                geom_atlas = atlas_fmri_geom if func_wt > 0 else  atlas_folding_geom
                if which_training == 'func':
                    input0 = mrisps_func_test[sno]
                    input1 = func_atlas.data[...,np.newaxis]
                    input2 = mrisps_geom_test[sno]
                else:
                    input0 = mrisps_geom_test[sno]
                    input1 = geom_atlas.data
                    input2 = mrisps_func_test[sno]

                # with strategy.scope():
                if 1:
                    model_inputs = [input0[np.newaxis], input1[np.newaxis], input2[np.newaxis]]
                    pred = model_semi.predict(model_inputs)
                
                    if which_training == 'func':
                        geom = pred[2]
                        func = pred[0]
                    else:
                        geom = pred[0]
                        func = pred[2]
                        
                    avg_geoms.append(geom[0,...])
                    avg_funcs.append(func[0,...,0])

            # end if/else for target/atlas

        # end subject for loop
        geom_list.append(avg_geoms)
        func_list.append(avg_funcs)
        mean_geom = np.array(avg_geoms).mean(axis=0).squeeze()
        std_geom = np.array(avg_geoms).std(axis=0).squeeze()
        mean_func = np.array(avg_funcs).mean(axis=0).squeeze()
        std_func = np.array(avg_funcs).std(axis=0).squeeze()

        # end of atlas/subject if
        # these are the summary maps for each different functional weighting
        mean_geoms.append(mean_geom)
        std_geoms.append(std_geom)
        mean_funcs.append(mean_func)
        std_funcs.append(std_func)

        if save_images:
            print('saving average images')
            fname = name + f'.func_wt_%2.2f' % 100
            fs.Image(mean_geom).write(fname + '.mean_geom.mrisp.mgz')
            fs.Image(std_geom).write(fname + '.std_geom.mrisp.mgz')
            fs.Image(mean_func).write(fname + '.mean_func.mrisp.mgz')
            fs.Image(std_func).write(fname + 'std_func.mrisp.mgz')

fv = fs.Freeview()
for fno, func_wt in enumerate(func_wt_list):
    sopts = ':visible=1:locked=1' if fno == 0 else ':visible=0:locked=1:linked=1'
    fopts = ':visible=1' if fno == len(func_wt_list)-1 else ':visible=0'
    fopts += ':heatscale=2,3:locked=0'
    fv.vol(mean_geoms[fno][pad:-pad, pad:-pad, np.newaxis, :], name=f'geom_fwt_{func_wt}', opts=sopts)
    fv.vol(mean_funcs[fno][pad:-pad, pad:-pad, ...], name=f'func_{func_wt}', colormap='heat',opts=fopts)

fv.show()


if build_atlases:
    for fno, func_wt in enumerate(func_wt_list):
        print('saving atlas images')
        fname = name + f'atlas.func_wt_%2.2f' % func_wt
        fs.Image(mean_geom).write(fname + '.mean_geom.mrisp.mgz')
        fs.Image(std_geom).write(fname + '.std_geom.mrisp.mgz')
        fs.Image(mean_func).write(fname + '.mean_func.mrisp.mgz')
        fs.Image(std_func).write(fname + 'std_func.mrisp.mgz')
    
t1 = 2
t2 = 3
fv = fs.Freeview()
geom_overlay = fsa_sphere.sample_parameterization(atlas_fmri_geom.data[pad:-pad, pad:-pad, 2])
func_overlay = fsa_sphere.sample_parameterization(atlas_fmri_func.data[pad:-pad, pad:-pad])
func_tag = fv.OverlayTag(func_overlay, name='func_atlas', threshold=f'{t1},{t2}')
geom_tag = fv.OverlayTag(geom_overlay, name='geom_atlas', threshold=f'{t1},{t2}')
fv.surf(fsa_inf, overlay=func_overlay, curvature=geom_overlay, opts=f":overlay_threshold={t1},{t2}:name=atlas:locked=1")
for ono, wt in enumerate(func_wt_list):
    geom_overlay = fsa_sphere.sample_parameterization(mean_geoms[ono][pad:-pad, pad:-pad, 2])
    func_overlay = fsa_sphere.sample_parameterization(mean_funcs[ono][pad:-pad, pad:-pad])
    func_tag = fv.OverlayTag(func_overlay, name='func_%2.2f' % wt, threshold=f'{t1},{t2}')
    geom_tag = fv.OverlayTag(geom_overlay, name='geom_%2.2f' % wt, threshold=f'{t1},{t2}')
    fv.surf(fsa_inf, overlay=func_overlay, curvature=geom_overlay, opts=f":overlay_threshold={t1},{t2}:name=func_{wt}")

#fv.surf(fsa_inf, overlay=func_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
#fv.surf(fsa_inf, overlay=func_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
fv.show(title='fs average langauge overlays', opts='--viewport 3D')
