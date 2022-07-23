import matplotlib.pyplot as plt
import pdb as gdb
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
import nibabel as nib
from sklearn.utils import class_weight
from nibabel import processing as nip
import numpy as np
import scipy.ndimage.morphology as morph
import surfa as sf
import os,socket
from netshape import *
from dipy.align.reslice import reslice
import neuron as ne
import voxelmorph as vxm
from netparms import *
from freesurfer import deeplearn as fsd
from freesurfer.deeplearn.utils import WeightsSaver, ModelSaver, utils, pprint
import imageio

cwd = os.getcwd()
host = socket.gethostname()
idir = '/autofs/cluster/lcnextdata1/ChestXray-NIHCC'
ndata = -1
gpu_number = 0
if (host == 'tebo.nmr.mgh.harvard.edu'):
    gpu_number = 0
elif (host == 'serena.nmr.mgh.harvard.edu'):
    gpu_number = 0
#    ndata = 5
elif (host == 'sulc.nmr.mgh.harvard.edu'):
    gpu_number = 1
elif (host == 'mlscgpu1'):
    gpu_number = 5
elif (host == 'mlscgpu2.nmr.mgh.harvard.edu'):
    gpu_number = 5

    
def BatchGenerator(ipath, images, batch_size=4,return_warp=False, bidir=True):
    nimages = len(images)
    im0 = imageio.imread(os.path.join(ipath, images[0]))
    warp = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([2]))
    batch_source = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    batch_target = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    found = 0
    while (True):
        ind1 = np.random.randint(0,nimages)
        ind2 = np.random.randint(0,nimages)
        im_src = imageio.imread(os.path.join(ipath, images[ind1]))
        im_trg = imageio.imread(os.path.join(ipath, images[ind2]))
        if len(im_src.shape) > 2:
            im_src = im_src[...,0]
        if len(im_trg.shape) > 2:
            im_trg = im_trg[...,0]
        batch_source[found,...,0] = im_src / im_src.max()
        batch_target[found,...,0] = im_trg / im_trg.max()
        found += 1
        if found >= batch_size:
            if return_warp == True:
                yield([batch_source, batch_target], [batch_target,warp])
            else:
                yield([batch_source, batch_target], batch_target)
            found = 0

train_affine = True
scale=1
print('running on host %s, GPU %d, train_affine %s' % (host, gpu_number, str(train_affine)))
fsd.configure(gpu=gpu_number)

fp = open(os.path.join(idir, 'train_val_list.txt'))
images = fp.read().split('\n')[2:-1];
fp.close()
nimages = len(images)
im = imageio.imread(os.path.join(idir, 'images', images[0]))

target_shape = (int(im.shape[0]/scale),int(im.shape[1]/scale))


batch_size=16

learning_rate = 0.001*.00001

rigid = False
which_loss = 'tukey'
which_loss = 'ncc'
which_loss = 'mse'
symmetrize = False
affine_blurs = [8,4,2,1]
affine_blurs = [1]
bidir = False
model = vxm.networks.VxmAffineDense(target_shape, enc_nf, dec_nf, enc_nf_affine=enc_nf_affine,rigid=rigid, affine_bidir=bidir, affine_blurs=affine_blurs)

ncc = vxm.losses.NCC([25,25])
tukey_c = .25
tukey = vxm.losses.TukeyBiweight(c=tukey_c)
if which_loss == 'ncc':
    aloss=ncc.loss
elif which_loss == 'mse':
    aloss=keras.losses.mse
elif which_loss == 'tukey':
    aloss = tukey.loss

ldir = 'logs/cxr.%s.train_affine.%s' % (which_loss, str(train_affine))
tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=ldir)

reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='loss', factor=0.5, patience=10, min_lr=1e-8, verbose=True, min_delta=0.0001)
affine_wt_fname = 'cxr.vxm.mpr.%s.%s.sym.%s.train_affine.%s.h5' % (str(rigid), which_loss, str(symmetrize),str(train_affine))
callbacks = [ModelSaver(model, 100, affine_wt_fname,cp_iters=10), reduce_lr, tensorboard_callback]
steps_per_epoch = int(128)

affine_model = model.references.affine_model


bg = BatchGenerator(os.path.join(idir, 'images'),images, batch_size=batch_size,return_warp=False, bidir=bidir)
affine_model.compile(optimizer=keras.optimizers.Adam(lr=.1*learning_rate), loss=aloss)
initial_epoch = 50
epochs = initial_epoch+50
#fithistr = affine_model.fit_generator(bg, steps_per_epoch = steps_per_epoch, epochs = epochs, verbose=1, callbacks=callbacks, class_weight=None, initial_epoch=initial_epoch)
if 0:
    np.random.seed(1)
    bg2 = BatchGenerator(os.path.join(idir, 'images'),images, batch_size=batch_size,return_warp=False, bidir=bidir)
    inb,outb = next(bg2)
    p = affine_model.predict(inb)
    fv = sf.vis.Freeview()
    fv.add_image(np.transpose(inb[0][...,0], (1,2,0)), name='src')
    fv.add_image(np.transpose(inb[1][...,0], (1,2,0)), name='trg')
    if isinstance(p, list):
        p = p[0]
    fv.add_image(np.transpose(p[...,0], (1,2,0)), name='p')
    fv.show()
    affine_pred_model = affine_model.get_predictor_model()
    affine = affine_pred_model.predict(inb)
    if isinstance(affine, list):
        affine = affine[0]
    maffine = np.reshape(np.append(affine[0,:],np.array([0,0,1])),(3,3))
    pprint(maffine)
    w,v = np.linalg.eig(maffine)
    print(w,v)
    assert 0

if train_affine == False:
    print('disabling training in affine layers')
    set_trainable(affine_model, False)

affine_model.save(affine_wt_fname)
        
ldir = 'logs/cxr.%s.train_affine.%s.nl' % (which_loss, str(train_affine))
tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=ldir)

nl_wt_fname = 'cxr.vxm.mpr.%s.%s.sym.%s.train_affine.%s.nl.h5' % (str(rigid), which_loss, str(symmetrize),str(train_affine))
which_loss = 'mse'
aloss=keras.losses.mse
losses = [aloss, vxm.losses.Grad().loss]
smoothness_weight=1
model.compile(optimizer=keras.optimizers.Adam(lr=.1*learning_rate), loss=losses, loss_weights=[1,smoothness_weight])

bg_nl = BatchGenerator(os.path.join(idir, 'images'),images, batch_size=batch_size,return_warp=True, bidir=False)

nl_epochs = 50
initial_epoch = 0
epochs = initial_epoch+nl_epochs
fithistr = model.fit_generator(bg_nl, steps_per_epoch = steps_per_epoch, epochs = initial_epoch+epochs, verbose=1, callbacks=callbacks, class_weight=None, initial_epoch=initial_epoch)
affine_model.save(affine_wt_fname)
model.save(nl_wt_fname)
if 1:
    np.random.seed(1)
    bg_nl2 = BatchGenerator(os.path.join(idir, 'images'),images, batch_size=batch_size,return_warp=True, bidir=False)

    inb,outb=next(bg_nl2)
    src = inb[0]
    trg = outb[0]
    nl_p, nl_warp = model.predict(inb)
    jlist = []
    affine_p = affine_model.predict(inb)
    if isinstance(affine_p, list):
        affine_p = affine_p[0]
    fv = sf.vis.Freeview()
    fv.add_image(np.transpose(src, (1,2,0,3)), name='src')
    fv.add_image(np.transpose(trg, (1,2,0,3)), name='trg')
    if which_loss == 'tukey':
        name = 'pred.rigid.%s.%s.sym.%s.affine.%s' % (str(rigid),str(which_loss)+str(tukey_c),str(symmetrize),str(train_affine))
    else:
        name = 'pred.rigid.%s.%s.sym.%s.affine.%s' % (str(rigid),str(which_loss),str(symmetrize),str(train_affine))
    fv.add_image(np.transpose(affine_p, (1,2,0,3)), name='affine'+name)
    fv.add_image(np.transpose(nl_p, (1,2,0,3)), name='NL'+name)
    fv.show()
    affine_pred_model = affine_model.get_predictor_model()
    affine = affine_pred_model.predict(inb)
    maffine = np.reshape(np.append(affine[0,:],np.array([0,0,1])),(3,3))
    pprint(maffine)
    w,v = np.linalg.eig(maffine)
    print(w,v)
