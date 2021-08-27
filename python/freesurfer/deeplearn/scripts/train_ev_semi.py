#%load_ext autoreload
#%autoreload 2

import os, socket, sys, copy
from glob import glob
from tensorflow import keras
from tensorflow.keras import backend as K
import tensorflow as tf

import freesurfer as fs
from freesurfer import deeplearn as fsd
import tensorflow as tf
import neurite_sandbox as nes
import voxelmorph as vxm

batch_size = 8

model_dir = 'models'
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

tf.config.get_visible_devices()

# paths to data and list of subjects and such
base_dir = '/autofs/cluster/p41/resting_state/Freesurfer_Surface'
sdir1 = 'CG1'
sdir2 = 'CG2'
fsdir = os.path.join(base_dir, sdir1, 'FS')
fsdir_val = os.path.join(base_dir, sdir2, 'FS')

hemi = 'lh'
sphere_name = '%s.sphere.rot' % hemi
sulc_name = '%s.sulc' % hemi

fsa_surf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.sphere.reg'))
fsa_inf = fs.Surface.read(os.path.join(fsdir, 'fsaverage', 'surf', hemi+'.inflated'))

contrast_name = 'SvsN'
func_dir = os.path.join('bold','bold.self.sm0.%s.lang' % hemi, contrast_name, 't.nii.gz')
spaths = glob(os.path.join(fsdir, 'sub???'))
crop = -1
spaths = spaths[0:crop]
spaths_val = glob(os.path.join(fsdir_val, 'sub???'))

# load all the data (if not previously loaded)
pad = 8
smooth_steps = 200  # apply Gaussian smoothing on the surface (sigma~sqrt(smoothing))
if 'spheres' not in locals() and 'spheres' not in globals():
    mrisps_geom, mrisps_func, spheres, snames = fsd.surf_utils.load_func_and_spheres(spaths, base_dir, func_dir, sdir1, sphere_name, hemi, pad=pad, smooth_steps=smooth_steps)
    mrisps_geom_val, mrisps_func_val, spheres_val, snames_val = fsd.surf_utils.load_func_and_spheres(spaths_val, base_dir, func_dir, sdir2, sphere_name, hemi, pad=pad, smooth_steps=smooth_steps)

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
fgen = fsd.surf_utils.mrisp_semi_gen(mrisps_func, mrisps_geom, batch_size=batch_size, use_rand=True, warp_downsize=warp_downsize)
ggen = fsd.surf_utils.mrisp_semi_gen(mrisps_geom, mrisps_func, batch_size=batch_size, use_rand=False, warp_downsize=warp_downsize)

which_training = 'sup'
which_training = 'semi'
which_semi = 'POW'
which_semi = 'L2'
which_semi = 'NCC'

if which_training == 'sup':
    gen = fgen  # function goes through convs
else:
    gen = ggen  # geometry goes through convs

ninputs = 1  # sulc
mrisp_shape = mrisps_geom[0].shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)
warp_sphere_loss = fsd.losses.spherical_loss(tuple(np.array(mrisp_shape_nopad)//warp_downsize), pad=pad//2)
sphere_loss = fsd.losses.spherical_loss(mrisp_shape_nopad, pad=pad, threshold=threshold, win=[41,41])
if which_semi == 'NCC':
    semi_loss = sphere_loss.NCC_loss(1)
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

warp_smooth_wt = .1
geom_match_wt = 1

fhists = []
mean_geoms = []
std_geoms = []
mean_funcs = []
std_funcs = []
geom_list = []
func_list = []

# vxm unet architecture
enc_nfeats = [64, 96, 128, 128]
dec_nfeats = [128, 64, 64, 64, 32,32]
unet_nfeats = [enc_nfeats, dec_nfeats]
lr=1e-4
func_wt_list = [0,.1,1,10,100]
for func_wt in func_wt_list:
    if func_wt > 1:
        func_wt_this_round = 1
        geom_match_wt_this_round = geom_match_wt / func_wt
    else:
        func_wt_this_round = func_wt
        geom_match_wt_this_round = 1

    print('training with %s loss weight %2.2f and which_semi = %s, geom_wt=%2.2f, warp_smooth_wt=%2.2f' % (which_training, func_wt_this_round,which_semi, geom_match_wt_this_round,warp_smooth_wt))

    glist = []
    for gno in range(2,4):
        glist.append('/gpu:%d' % gno)
        
    strategy = tf.distribute.MirroredStrategy()
    
    with strategy.scope():
        model_semi = vxm.networks.VxmDenseSemiSupervisedSeg(
            mrisp_shape, 
            1, 
            nb_unet_features=unet_nfeats, 
            int_downsize=warp_downsize, 
            seg_downsize=1)

        # callbacks
        model_fname = '%s/%s.s2.surf2surf.%s%2.2f.warp_smooth_wt%2.2f.%s.h5' % \
                      (model_dir, hemi,which_training,func_wt,warp_smooth_wt,which_semi)
        stopping_callback = tf.keras.callbacks.EarlyStopping(
            monitor='loss', mode='min', verbose=1, patience=200)
        mc = tf.keras.callbacks.ModelCheckpoint(
            model_fname, monitor='loss', mode='min', save_best_only=True)
        lr_callback = tf.keras.callbacks.ReduceLROnPlateau(
            factor=0.8, patience=100, cooldown=10, monitor='loss')
        callbacks = [mc, lr_callback, stopping_callback]

        # losses and compile
        losses = [semi_loss, warp_sphere_loss.gradientLoss(penalty='l2'), semi_loss]
        if which_training == 'sup':
            loss_weights = [func_wt_this_round,warp_smooth_wt,geom_match_wt_this_round]
        else:
            loss_weights = [geom_match_wt_this_round,warp_smooth_wt,func_wt_this_round]

        nes.utils.check_and_compile(model_semi, gen, optimizer=keras.optimizers.Adam(lr=lr), 
                                    loss=losses, loss_weights=loss_weights)
        fhist = model_semi.fit(gen, epochs=2000, steps_per_epoch=100, callbacks=callbacks)

    fhists.append(fhist)

    plt.figure()
    plt.plot(fhist.history['vxm_dense_transformer_loss'])
    plt.plot(fhist.history['seg_transformer_loss'])
    plt.plot(fhist.history['loss'])
    plt.grid()
    if which_training == 'sup':
        plt.title('func = direct, geom = semi (wt=%2.2f)' % func_wt)
        plt.legend(['func', 'geom', 'loss'])
    else:
        plt.title('geom = direct, func = semi (wt=%2.2f)' % func_wt)
        plt.legend(['geom', 'func', 'loss'])

    plt.show(block=False)

    if 0:
        model_fsemi = vxm.networks.VxmDenseSemiSupervisedSeg(mrisp_shape, 1, nb_unet_features=unet_nfeats, int_downsize=1, seg_downsize=1)
        model_fsemi.compile(optimizer=keras.optimizers.Adam(lr=1.0*lr), loss=losses, loss_weights=loss_weights)
        fhist2 = model_fsemi.fit(fgen, epochs=200, steps_per_epoch=100, callbacks=callbacks)
        plt.figure()
        plt.plot(fhist2.history['transformer_loss'])
        plt.plot(fhist2.history['seg_transformer_loss'])
        plt.grid()
        plt.title('func_wt = %2.2f' % (func_wt))
        plt.title('func = direct, geom = semi')
        plt.legend(['func', 'geom'])
        plt.show(block=False)

    model_semi.save_weights(model_fname)

    if 0:
        inb, outb = next(gen)
        pred = model_semi.predict(inb)
        moving = np.transpose(inb[0], (1,2,3,0))
        fixed = np.transpose(inb[1], (1,2,3,0))
        moving_func = np.transpose(inb[2], (1,2,3,0))
        fixed_func = np.transpose(outb[2], (1,2,3,0))
        moving_warped = np.transpose(pred[0], (1,2,3,0))
        moving_func_warped = np.transpose(pred[2], (1,2,3,0))
        fv = fs.Freeview()
        fv.vol(fixed, name='fixed')
        fv.vol(moving, name='moving', visible='0')
        fv.vol(moving_warped, name='moving_warped')
        fv.vol(moving_func, name='func_moving', colormap='heat', visible='0')
        fv.vol(moving_func_warped, name='moving_func_warped_%2.2f' % func_wt, colormap='heat', visible='0')
        fv.vol(fixed_func, name='func_fixed', colormap='heat', visible='0')
        fv.vol((moving_func+fixed_func)/2, name='avg_b4', colormap='heat', visible='1')
        fv.vol((moving_func_warped+fixed_func)/2, name='avg_after', colormap='heat', visible='1')
        fv.show(title='func_wt=%2.1f' % func_wt)

    warped_funcs = []
    warped_geoms = []
    avg_warps = []
    avg_geoms = []
    avg_funcs = []
    trans_model = vxm.networks.Transform(mrisp_shape)
    registration_model = model_semi.get_registration_model()
    print('computing subject warp matrix')
    for sno, sname in enumerate(tqdm(snames)):
        wgeoms = []
        inputs = [[],[],[]]
        for sno2, sname2 in enumerate(snames):
#            inputs = [mrisps_geom[sno][np.newaxis], mrisps_geom[sno2][np.newaxis], mrisps_func[sno][np.newaxis]]
            if which_training == 'sup':
                inputs[0].append(mrisps_func[sno])
                inputs[1].append(mrisps_func[sno2])
                inputs[2].append(mrisps_geom[sno])
            else:
                inputs[0].append(mrisps_geom[sno])
                inputs[1].append(mrisps_geom[sno2])
                inputs[2].append(mrisps_func[sno])

        model_inputs = [np.array(inputs[0]), np.array(inputs[1]), np.array(inputs[2])]
        pred = model_semi.predict(model_inputs)
        wgeoms = pred[0]
        warps = pred[1]
        wfuncs = pred[2]
        warps = registration_model.predict(model_inputs)

        warped_funcs.append(wfuncs)
        warped_geoms.append(wgeoms)
                         
        avg_warp = np.array(warps).mean(axis=0).squeeze()
        avg_geom = np.array(wgeoms).mean(axis=0).squeeze()
        avg_func = np.array(wfuncs).mean(axis=0).squeeze()
        avg_warps.append(avg_warp)
        warped_geom = trans_model.predict([mrisps_geom[sno][np.newaxis], avg_warp[np.newaxis]])
        warped_func = trans_model.predict([mrisps_func[sno][np.newaxis], avg_warp[np.newaxis]])

        avg_geoms.append(warped_geom)
        avg_funcs.append(warped_func)
  
    fv = fs.Freeview()
    fv.vol(np.transpose(np.array(avg_geoms[0:20]), (2,3,4,0,1)), name='geom')
    fv.vol(np.transpose(np.array(avg_funcs[0:20]), (2,3,4,0,1)), name='func_%2.2f' % func_wt, colormap='heat',opts=':heatscale=2,3')
    fv.show()

    geom_list.append(avg_geoms)
    func_list.append(avg_funcs)

    mean_geom = np.array(avg_geoms).mean(axis=0).squeeze()
    std_geom = np.array(avg_geoms).std(axis=0).squeeze()
    mean_func = np.array(avg_funcs).mean(axis=0).squeeze()
    std_func = np.array(avg_funcs).std(axis=0).squeeze()

    mean_geoms.append(mean_geom)
    std_geoms.append(std_geom)
    mean_funcs.append(mean_func)
    std_funcs.append(std_func)

    fs.Image(mean_geom).write('%s.d%d.func_wt%2.2f.warp_wt%2.2f.smooth%d.mean_geom.mrisp.mgz' % (hemi, warp_downsize, func_wt, warp_smooth_wt, smooth_steps))
    fs.Image(std_geom).write( '%s.d%d.func_wt%2.2f.warp_wt%2.2f.smooth%d.std_geom.mrisp.mgz' % (hemi, warp_downsize, func_wt, warp_smooth_wt, smooth_steps))
    fs.Image(mean_func).write('%s.d%d.func_wt%2.2f.warp_wt%2.2f.smooth%d.mean_func.mrisp.mgz' % (hemi, warp_downsize,func_wt, warp_smooth_wt, smooth_steps))
    fs.Image(std_func).write('%s.d%d.func_wt%2.2f.warp_wt%2.2f.smooth%d.std_func.mrisp.mgz' % (hemi, warp_downsize, func_wt, warp_smooth_wt, smooth_steps))


geom_tag_list = []
func_tag_list = []
fv = fs.Freeview()
for ono, wt in enumerate(func_wt_list):
    geom_overlay = fsa_surf.sample_parameterization(mean_geoms[ono][pad:-pad,:])
    func_overlay = fsa_surf.sample_parameterization(mean_funcs[ono][pad:-pad,:])
    func_tag = fv.OverlayTag(func_overlay, name='func_%2.2f' % wt, threshold='2,3')
    geom_tag = fv.OverlayTag(geom_overlay, name='geom_%2.2f' % wt, threshold='2,3')
    func_tag_list.append(func_tag)
    geom_tag_list.append(geom_tag)


#fv.surf(fsa_inf, overlay=func_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
fv.surf(fsa_inf, overlay=geom_tag_list, curvature=geom_overlay, opts=":overlay_threshold=2,3")
fv.show()

