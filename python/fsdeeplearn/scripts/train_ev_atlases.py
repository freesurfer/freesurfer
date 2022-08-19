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
import fsdeeplearn as fsd
import tensorflow as tf
import neurite_sandbox as nes
import voxelmorph as vxm

from neurite_sandbox.tf.utils.utils import plot_fit_callback as pfc

save_images = False
save_images = True

dofit = True
dofit = False

which_training = 'geom'
which_training = 'func'

which_semi = 'POW'
which_semi = 'L2'
which_semi = 'NCC'

batch_size = 8

warp_smooth_wt = .01
warp_smooth_wt = .1
geom_match_wt = .1
geom_match_wt = 1

hemi = 'lh'

smooth_steps = 200  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))
smooth_steps = 50  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))
smooth_steps = 1  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))

model_dir = 'models'
name = f'{model_dir}/{hemi}.atlas.{which_training}.warp_smooth_wt_%2.2f.geom_wt_{geom_match_wt}.smooth_{smooth_steps}.{which_semi}' % (warp_smooth_wt)

gpu_id = -1
host = socket.gethostname()
if os.getenv('NGPUS'):
    ngpus = int(os.getenv('NGPUS'))
#    ngpus=1
    gpu_str = '0'
    for g in range(1,ngpus):
        gpu_str += ',%d' % g
    os.environ["CUDA_VISIBLE_DEVICES"] = gpu_str
    print('reading %d GPUS from env and setting CUDA_VISIBLE_DEVICES to %s' % (ngpus, gpu_str))
elif os.getenv('CUDA_VISIBLE_DEVICES'):
    gpu_list = os.getenv('CUDA_VISIBLE_DEVICES')
    ngpus = len(gpu_list.split(','))
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

print(f'dofit {dofit}, which_semi {which_semi}, which_training {which_training}, smooth_steps {smooth_steps}')

tf.config.get_visible_devices()

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

if 'spheres' not in locals() and 'spheres' not in globals() and smooth_steps != smooth_read:
    geom_fname = f'mrisps_geom_smooth_{smooth_steps}.npy'
    func_fname = f'mrisps_func_smooth_{smooth_steps}.npy'
    spheres_fname = f'spheres_smooth_{smooth_steps}.npy'
    snames_fname = f'snames_smooth_{smooth_steps}.npy'
    geom_fname_val = f'mrisps_geom_smooth_{smooth_steps}_val.npy'
    func_fname_val = f'mrisps_func_smooth_{smooth_steps}_val.npy'
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

round = 0
if round == 0:
    geom_atlas = np.array(mrisps_geom).mean(axis=0).squeeze()
    func_atlas = np.array(mrisps_func).mean(axis=0).squeeze()
    fgen = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_func, mrisps_geom, func_atlas[...,np.newaxis], 
                                               geom_atlas, batch_size=batch_size, 
                                               use_rand=True, warp_downsize=warp_downsize)
    ggen = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, geom_atlas, 
                                               func_atlas[...,np.newaxis], batch_size=batch_size, 
                                               use_rand=True, warp_downsize=warp_downsize)
else:
    fname = 'models/lh.surf2surf.sup.warp_smooth_wt_0.10.smooth_1.NCC.func_wt_100.00'

    semi_name = 'models/lh.surf2surf.semi.warp_smooth_wt_0.10.geom_wt_1.smooth_1.NCC.func_wt_0.00.h5'
    sup_name = 'models/lh.surf2surf.sup.warp_smooth_wt_0.10.geom_wt_1.smooth_1.NCC.func_wt_10.00.h5'
    print('saving atlas images')
    fname_func = aname + f'atlas.func_wt_%2.2f' % 0
    fname_geom = aname + f'atlas.func_wt_%2.2f' % 0

    print(f'loading atlases from {fname_geom} and {fname_func}')
    func_atlas_func = fs.Image.read(fname_func + '.mean_func.mrisp.mgz')
    func_atlas_geom = fs.Image.read(fname_func + '.mean_geom.mrisp.mgz')
    geom_atlas_func = fs.Image.read(fname_geom + '.mean_func.mrisp.mgz')
    geom_atlas_geom = fs.Image.read(fname_geom + '.mean_geom.mrisp.mgz')

    # use atlases built from geometry-only
    fgen0 = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_func, mrisps_geom, geom_atlas_func.data[...,np.newaxis], 
                                               geom_atlas_geom.data, batch_size=batch_size, 
                                               use_rand=True, warp_downsize=warp_downsize)
    ggen0 = fsd.surf_utils.mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, geom_atlas_geom.data, 
                                                geom_atlas_func.data[...,np.newaxis], batch_size=batch_size, 
                                                use_rand=True, warp_downsize=warp_downsize)

if which_training == 'geom':
    ninputs = mrisps_geom[0].shape[-1] 
    semi_inputs = 1
else:
    semi_inputs = mrisps_geom[0].shape[-1] 
    ninputs = 1
    
mrisp_shape = mrisps_geom[0].shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)
warp_sphere_loss = fsd.losses.spherical_loss(tuple(np.array(mrisp_shape_nopad)//warp_downsize), pad=pad//2)
sphere_loss = fsd.losses.spherical_loss(mrisp_shape_nopad, pad=pad, threshold=threshold, win=[21,21])
if which_semi == 'NCC':
    semi_loss = sphere_loss.NCC_signed_loss(weight=1)
elif which_semi == 'L2':
    semi_loss = sphere_loss.l2_loss(1, None)
else:
    semi_loss = sphere_loss.power_loss

if 0:
    inb,outb = next(gen)
    pred = model_semi.predict(inb)
    t1 = tf.convert_to_tensor(inb[1], dtype=tf.float32)
    t2 = tf.convert_to_tensor(pred[0], dtype=tf.float32)
    t3 = tf.convert_to_tensor(outb[0], dtype=tf.float32)
    l = semi_loss(t1, t2)
    wl = warp_sphere_loss.gradientLoss(penalty='l2')(_, pred[1])
    print(K.eval(l))
    assert 0

fhists = []

# vxm unet architecture
enc_nfeats = [64, 96, 128, 128]
dec_nfeats = [128, 64, 64, 64, 32,32]
unet_nfeats = [enc_nfeats, dec_nfeats]
lr = 1e-4
func_wt_list = [0,.1,1,10,100]
func_wt_list = [.1]
func_wt_list = [100]
func_wt_list = [1]
func_wt_list = [10]
func_wt_list = [0]
func_wt_list = [0, 10]
if not dofit:
    if which_training == 'func':
        func_wt_list = [100]
    else:
        func_wt_list = [0, 0.1, 1, 10, 100]
        func_wt_list = [0, 1, 10, 100]

func_wt = 0 if which_training == 'geom' else 10

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
#with strategy.scope():
if 1:
    model_vxm = vxm.networks.VxmDenseSemiSupervisedSeg(
        mrisp_shape, 
        semi_inputs, 
        src_feats=ninputs,
        trg_feats=ninputs,
        nb_unet_features=unet_nfeats, 
        int_resolution=warp_downsize, 
        seg_resolution=1)
        
    if dofit:
        print('training with func loss weight %2.2f and which_semi = %s, geom_wt=%2.2f, warp_smooth_wt=%2.2f' % 
              (func_wt_this_round,which_semi, geom_match_wt_this_round,warp_smooth_wt))

        new_out1 = KL.Lambda(lambda x : x, name='warp')(model_vxm.outputs[1])
        if which_training == 'func':
            new_out0 = KL.Lambda(lambda x : x, name='fmri')(model_vxm.outputs[0])
            new_out2 = KL.Lambda(lambda x : x, name='geom')(model_vxm.outputs[2])
        else:
            new_out0 = KL.Lambda(lambda x : x, name='geom')(model_vxm.outputs[0])
            new_out2 = KL.Lambda(lambda x : x, name='fmri')(model_vxm.outputs[2])
            
        model_semi = tf.keras.Model(model_vxm.inputs, [new_out0, new_out1, new_out2])
        model_semi.references = model_vxm.references
        # losses and compile
        losses = [ semi_loss, warp_sphere_loss.gradientLoss(penalty='l2'), semi_loss]
        if which_training == 'func':
            loss_weights = [func_wt_this_round,warp_smooth_wt,geom_match_wt_this_round]
            gen = fgen
        else:   # training only geometry through convs
            loss_weights = [geom_match_wt_this_round,warp_smooth_wt,func_wt_this_round]
            gen = ggen

        nes.utils.check_and_compile(model_semi, gen, optimizer=keras.optimizers.Adam(lr=lr), 
                                    loss=losses, loss_weights=loss_weights)
        print(f'writing fit history to {write_cb.fname} and model weights to {lr_cb.fname}')
        callbacks = [write_cb, lr_cb]
        fhist = model_semi.fit(gen, epochs=10000, steps_per_epoch=200, callbacks=callbacks)
        # fhist = model_semi.fit(gen, epochs=10000, steps_per_epoch=200, callbacks=callbacks, initial_epoch=5000)
        
        fhists.append(fhist)
        # model_semi.save_weights(model_fname)
    else:
        vxm_geom = vxm.networks.VxmDenseSemiSupervisedSeg(
            mrisp_shape, 
            1, 
            src_feats=mrisps_geom[0].shape[-1],
            trg_feats=mrisps_geom[0].shape[-1],
            nb_unet_features=unet_nfeats, 
            int_resolution=warp_downsize, 
            seg_resolution=1)
        vxm_func = vxm.networks.VxmDenseSemiSupervisedSeg(
            mrisp_shape, 
            mrisps_geom[0].shape[-1], 
            src_feats=1,
            trg_feats=1,
            nb_unet_features=unet_nfeats, 
            int_resolution=warp_downsize, 
            seg_resolution=1)


func_wt_list = [0, 10]
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
    fit_names = []
    for fno, func_wt in enumerate(func_wt_list):    # for every different model

        if func_wt > 0:
            fname = f'{model_dir}/{hemi}.atlas.func.warp_smooth_wt_%2.2f.geom_wt_{geom_match_wt}.smooth_{smooth_steps}.{which_semi}' % (warp_smooth_wt)
            model_semi = vxm_func if func_wt > 0 else vxm_geom
            model_name = fname + f'.func_wt_%2.2f.h5' % func_wt
        else:
            fname = f'{model_dir}/{hemi}.atlas.geom.warp_smooth_wt_%2.2f.geom_wt_{geom_match_wt}.smooth_{smooth_steps}.{which_semi}' % (warp_smooth_wt)
            model_semi = vxm_geom
            model_name = fname + f'.func_wt_%2.2f.h5' % func_wt
    
        model_semi.load_weights(model_name)
        fit_names.append(".".join(model_name.split('.')[0:-1])+'.txt')
        warped_geoms = []
        warped_funcs = []
        trans_geom_model = vxm.networks.Transform(mrisp_shape, nb_feats=mrisps_geom[0].shape[-1])
        trans_func_model = vxm.networks.Transform(mrisp_shape)
        registration_model = model_semi.get_registration_model()
        print(f'warping each subject to the atlas for func_wt {func_wt}: {fno+1} of {len(func_wt_list)}')
        for sno, sname in enumerate(tqdm(snames_test)):   # compute and apply warp to implicit midspace for each subject
            if func_wt > 0:
                input0 = mrisps_func_test[sno]
                input1 = func_atlas[...,np.newaxis]
                input2 = mrisps_geom_test[sno]
            else:
                input0 = mrisps_geom_test[sno]
                input1 = geom_atlas
                input2 = mrisps_func_test[sno]

            # with strategy.scope():
            if 1:
                model_inputs = [input0[np.newaxis], input1[np.newaxis], input2[np.newaxis]]
                pred = model_semi.predict(model_inputs)
                
                if func_wt > 0:
                    geom = pred[2]
                    func = pred[0]
                else:
                    geom = pred[0]
                    func = pred[2]
                        
                warped_geoms.append(geom[0,...])
                warped_funcs.append(func[0,...,0])

        # end subject for loop
        geom_list.append(warped_geoms)
        func_list.append(warped_funcs)
        mean_geom = np.array(warped_geoms).mean(axis=0).squeeze()
        std_geom = np.array(warped_geoms).std(axis=0).squeeze()
        mean_func = np.array(warped_funcs).mean(axis=0).squeeze()
        std_func = np.array(warped_funcs).std(axis=0).squeeze()

        # end of atlas/subject if
        # these are the summary maps for each different functional weighting
        mean_geoms.append(mean_geom)
        std_geoms.append(std_geom)
        mean_funcs.append(mean_func)
        std_funcs.append(std_func)

    # end func_wt loop


if save_images:
    print('saving average images')
    fname = hemi + f'.{which_semi}.atlas.folding'
    fs.Image(mean_geoms[0]).write(fname + '.mean_geom.mrisp.mgz')
    fs.Image(std_geoms[0]).write(fname + '.std_geom.mrisp.mgz')
    fs.Image(mean_funcs[0]).write(fname + '.mean_func.mrisp.mgz')
    fs.Image(std_funcs[0]).write(fname + '.std_func.mrisp.mgz')

    fname = hemi + f'.{which_semi}.atlas.fmri'
    fs.Image(mean_geoms[1]).write(fname + '.mean_geom.mrisp.mgz')
    fs.Image(std_geoms[1]).write(fname + '.std_geom.mrisp.mgz')
    fs.Image(mean_funcs[1]).write(fname + '.mean_func.mrisp.mgz')
    fs.Image(std_funcs[1]).write(fname + '.std_func.mrisp.mgz')

fv = fs.Freeview()
for fno, func_wt in enumerate(func_wt_list):
    sopts = ':visible=1:locked=1' if fno == 0 else ':visible=0:locked=1'
    fopts = ':visible=1' if fno == len(func_wt_list)-1 else ':visible=0'
    fopts += ':heatscale=2,3:locked=0'
    fv.vol(mean_geoms[fno][pad:-pad, pad:-pad], name=f'geom_fwt_{func_wt}', opts=sopts)
    fv.vol(mean_funcs[fno][pad:-pad, pad:-pad], name=f'func_{func_wt}', colormap='heat',opts=fopts)

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
geom_overlay = fsa_sphere.sample_parameterization(geom_atlas[pad:-pad, pad:-pad,-1])
func_overlay = fsa_sphere.sample_parameterization(func_atlas[pad:-pad, pad:-pad])
func_tag = fv.OverlayTag(func_overlay, name='func_atlas', threshold=f'{t1},{t2}')
geom_tag = fv.OverlayTag(geom_overlay, name='geom_atlas', threshold=f'{t1},{t2}')
fv.surf(fsa_inf, overlay=func_overlay, curvature=geom_overlay, opts=f":overlay_threshold={t1},{t2}:name=atlas")
for ono, wt in enumerate(func_wt_list):
    geom_overlay = fsa_sphere.sample_parameterization(mean_geoms[ono][pad:-pad, pad:-pad,-1])
    func_overlay = fsa_sphere.sample_parameterization(mean_funcs[ono][pad:-pad, pad:-pad])
    func_tag = fv.OverlayTag(func_overlay, name='func_%2.2f' % wt, threshold=f'{t1},{t2}')
    geom_tag = fv.OverlayTag(geom_overlay, name='geom_%2.2f' % wt, threshold=f'{t1},{t2}')
    fv.surf(fsa_inf, overlay=func_overlay, curvature=geom_overlay, opts=f":overlay_threshold={t1},{t2}:name=func_{wt}")

#fv.surf(fsa_inf, overlay=func_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
#fv.surf(fsa_inf, overlay=func_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
fv.show(title='fs average langauge overlays', opts='--viewport 3D')
