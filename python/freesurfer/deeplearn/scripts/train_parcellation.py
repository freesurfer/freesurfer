
# system imports
import os, socket, glob, copy
from tensorflow import keras
import tensorflow as tf
from tensorflow.keras import backend as K
from tqdm import tqdm

# local imports
import neurite_sandbox as nes
import neurite as ne
import voxelmorph as vxm

# freesurfer imports
import freesurfer as fs
from freesurfer import deeplearn as fsd
from freesurfer.deeplearn.surf_utils import loadSphere
from freesurfer.deeplearn.utils import rebase_labels

pad = 16

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
elif os.getenv('CUDA_VISIBLE_DEVICES'):  # need to implement mirrored strategy for this
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



batch_size = 8
print('processing %s' % hemi)
sphere_name='sphere.rot'
inflated_name = 'inflated'
curv_name='sulc'
adir = os.path.join(os.getenv('FREESURFER_HOME'), 'average')
fs_atlas_fname = os.path.join(adir, hemi+'.average.curvature.filled.buckner40.tif')
tf.config.get_visible_devices()
feats = ['inflated.H', 'sulc', 'curv']


dkt_dir = os.path.join(atlas_dir, 'DKTatlas')
parc_name = 'labels.DKT31.manual.annot'
parc_name = 'labels.DKT31.manual.2.annot'
parc_fnames = glob.glob(os.path.join(dkt_dir, '*/label', hemi+'.'+parc_name))
ntraining = 50
nval = 10
ntest = len(parc_fnames) - (ntraining+nval)

# read in the data
# mrisp_geom holds the geometric features in spherical parameterization
# mrisp_annot are the maps of indices
# mrisp_onehot are the annots one-hot encoded
if 'read_hemi' not in locals() and 'read_hemi' not in globals():
    read_hemi = hemi

# check to see if data was loaded, and if not load it in
if read_hemi is not hemi or ('parcs' not in locals() and 'parcs' not in globals()):
    # load in all the surfaces and parcellations for training data
    dkt_subjects = [parc_fname.split('/')[-3] for parc_fname in parc_fnames]
    dkt_sphere_fnames = [os.path.join(dkt_dir, subject, 'surf', hemi+'.'+sphere_name) for subject in dkt_subjects]
    dkt_inflated_fnames = [os.path.join(dkt_dir, subject, 'surf', hemi+'.'+inflated_name) for subject in dkt_subjects]
    dkt_spheres = [fs.Surface.read(fname) for fname in tqdm(dkt_sphere_fnames)]
    dkt_inflated = [fs.Surface.read(fname) for fname in tqdm(dkt_inflated_fnames)]

    # build lookup table to compact labels into consecutive indices
    parcs = [fs.Overlay.read(parc_fname) for parc_fname in tqdm(parc_fnames)]
    all_labels = nes.py.utils.flatten_lists([list(np.where(parc.data<0, 0, parc.data)) for parc in parcs])
    al = [l if l >= 0 else 0 for l in all_labels]
    lab_to_ind, ind_to_lab = rebase_labels(al) # compact the label list
    nlabels = len(ind_to_lab)    # number of unique labels after compacting
    for parc in parcs:
        parc.data = lab_to_ind[np.where(parc.data<0, 0, parc.data)]

    # read in and spherically parameterize the geometry and parcellations
    dkt_mrisps_geom = []
    dkt_mrisps_annot = []
    for sno, subject in enumerate(tqdm(dkt_subjects)):
        for fno, curv_name in enumerate(feats): # geometry has 3 scalar fields
            surf_fname = os.path.join(dkt_dir, subject, 'surf', hemi+'.'+sphere_name)
            curv_fname = os.path.join(dkt_dir, subject, 'surf', hemi+'.'+curv_name)
            mrisp_tmp = loadSphere(surf_fname, curv_fname, padSize=pad)
            if fno == 0:  # first overlay - create mrisp
                mrisp = mrisp_tmp[...,np.newaxis]
            else:    # every other one concat to the end
                mrisp = np.concatenate((mrisp, mrisp_tmp[...,np.newaxis]), axis=-1)
        dkt_mrisps_geom.append(mrisp)
        mrisp_annot = np.transpose(dkt_spheres[sno].parameterize(parcs[sno], interp='nearest'),(1,0))
        if pad > 0:
            mrisp_annot = np.pad(mrisp_annot, ((pad,pad),(0,0)),'wrap')
            mrisp_annot = np.pad(mrisp_annot, ((0,0),(pad,pad)),'reflect')

        dkt_mrisps_annot.append(mrisp_annot)

    dkt_mrisps_onehot = [fsd.utils.np_one_hot(annot,nlabels) for annot in dkt_mrisps_annot]
    # split the data into training, validation and test
    dkt_mrisps_annot_training = dkt_mrisps_annot[0:ntraining]
    dkt_mrisps_onehot_training = dkt_mrisps_onehot[0:ntraining]
    dkt_mrisps_geom_training = dkt_mrisps_geom[0:ntraining]

    dkt_mrisps_annot_val = dkt_mrisps_annot[ntraining:ntraining+nval]
    dkt_mrisps_onehot_val = dkt_mrisps_onehot[ntraining:ntraining+nval]
    dkt_mrisps_geom_val = dkt_mrisps_geom[ntraining:ntraining+nval]


    dkt_mrisps_annot_test = dkt_mrisps_annot[ntraining+nval:]
    dkt_mrisps_onehot_test = dkt_mrisps_onehot[ntraining+nval:]
    dkt_mrisps_geom_test = dkt_mrisps_geom[ntraining+nval:]
    print('using %d labels' % nlabels)


mrisp_shape = dkt_mrisps_geom[0].shape[0:2]
mrisp_shape_nopad = (mrisp_shape[0]-2*pad, mrisp_shape[1]-2*pad)


# create the unet
nb_unet_features = [[64,64, 64], [64,64], [64,64], [64,64], [64,64]]

nfeats=64
conv_per_level=1
nb_unet_features = [[nfeats, nfeats, nfeats, nfeats, nfeats, nfeats, 2*nfeats, 2*nfeats, 2*nfeats, 2*nfeats], [2*nfeats, 2*nfeats, 2*nfeats, 2*nfeats, nfeats, nfeats, nfeats,nfeats, nfeats, nfeats,nfeats]]

unet = ne.models.unet(nb_unet_features, mrisp_shape+(len(feats),), None, 3, len(ind_to_lab), feat_mult=None, nb_conv_per_level=nb_conv_per_level)

# training, validation and testing generators
val_gen = fsd.surf_utils.parc_gen(dkt_mrisps_geom_val, dkt_mrisps_onehot_val, batch_size=1)
test_gen = fsd.surf_utils.parc_gen(dkt_mrisps_geom_test, dkt_mrisps_onehot_test, batch_size=ntest)
train_gen = fsd.surf_utils.parc_gen(dkt_mrisps_geom_training, dkt_mrisps_onehot_training, batch_size=8)

# loss and compilation
lr = 5e-4
sphere_loss = fsd.losses.spherical_loss(mrisp_shape_nopad[0:2], pad=pad, threshold=0)
dice_loss = sphere_loss.dice_loss(1)
optimizer = tf.keras.optimizers.Adam(lr=lr)

nes.utils.check_and_compile(unet,gen=train_gen,optimizer=optimizer, loss=[dice_loss])

# callbacks and training
reducelr = keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.75, patience=20, cooldown=1, min_lr=1e-7, verbose=1)
printlr = ne.callbacks.LRLog()
fit_results_fname = 'parc_unet.txt'
write_cb = nes.callbacks.WriteHist(fit_results_fname)
callbacks = [reducelr, write_cb, printlr]

hist = unet.fit(train_gen, epochs=250, steps_per_epoch=100, callbacks=callbacks, validation_data=val_gen, validation_steps=1)
unet.save(f'{hemi}.aparc.h5')


# plot training curves. This can be done during training due to the write_cb
nes.utils.plot_fit_callback(fit_results_fname, keys=['loss', 'val_loss','lr'], 
                            close_all=True, smooth=3)


# show results
inb, outb = next(test_gen)
pred = unet.predict(inb)

# show the spherical maps in freeview
fv = fs.Freeview(swap_batch_dim=True)
fv.vol(inb[...,1].squeeze(), name='inputs', opts=':linked=1:locked=1')
fv.vol(np.argmax(outb, axis=-1), name='targets', colormap='lut', opts=':linked=1:visible=0')
fv.vol(np.argmax(pred, axis=-1), name='pred', colormap='lut', opts=':linked=1,visible=1')
fv.show('unet_parcellation_results', verbose=True)

# now show them on the individual surfaces
aparcs = []
mparcs = []
sulcs = []
for sno, sphere in enumerate(tqdm(dkt_spheres[ntraining+nval:])):
    mrisp = np.argmax(pred[sno,pad:-pad,pad:-pad,...], axis=-1)
    aparc = sphere.sample_parameterization(ind_to_lab[np.transpose(mrisp)],interp='nearest')
    aparcs.append(aparc)

    mrisp = np.argmax(outb[sno,pad:-pad,pad:-pad,...], axis=-1)
    mparc = sphere.sample_parameterization(ind_to_lab[np.transpose(mrisp)],interp='nearest')
    mparcs.append(mparc)

    mrisp = dkt_mrisps_geom[sno][pad:-pad,pad:-pad,1]
    sulc = sphere.sample_parameterization(np.transpose(mrisp))
    sulcs.append(sulc)

num_to_show = 1
fv = fs.Freeview()
for n in range(num_to_show):
    ind = ntraining+nval+n
    overlay = fs.Overlay(aparcs[n])
    overlay.lut = parcs[0].lut
    mparc = copy.copy(parcs[ind])
    mparc.data = ind_to_lab[mparc.data]
    fv.surf(dkt_inflated[ind], annot=[overlay, mparc], curvature=sulcs[n], opts=':visible=%d' % (n==num_to_show-1))

# Ctrl + Alt + A will toggle between parcellations
fv.show(title='surface_parcellations', verbose=True)
