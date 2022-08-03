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
from fsdeeplearn.utils import WeightsSaver, ModelSaver, utils, pprint
import imageio, pydicom, gdcm, load_serial_cxr
from load_serial_cxr import load_timepoints
from neuron.utils import resize


bdir = '/autofs/cluster/lcnextdata1/CCDS_CXR/CXR-Serial/def_20200414'



cwd = os.getcwd()
host = socket.gethostname()
idir = '/autofs/cluster/lcnextdata1/ChestXray-NIHCC'
ndata = -1
gpu_number = 0
if (host == 'tebo.nmr.mgh.harvard.edu'):
    gpu_number = 1
elif (host == 'serena.nmr.mgh.harvard.edu'):
    gpu_number = 0
#    ndata = 5
elif (host == 'sulc.nmr.mgh.harvard.edu'):
    gpu_number = 1
elif (host == 'mlscgpu1'):
    gpu_number = 5
elif (host == 'mlscgpu2.nmr.mgh.harvard.edu'):
    gpu_number = 5


target_shape = (1024,1024)
vol_list, seg_list, dtrans_list, dcm_list, study_list, subject_name_list = load_timepoints(bdir, target_shape, ndilations=1)

print('%d training subjects found' % len(vol_list))
def BatchGenerator(vol_list, dtrans_list, batch_size=4,return_warp=False, bidir=False):
    nsubjects = len(vol_list)
    im0 = vol_list[0][0]
    d0 = dtrans_list[0][0]
    warp = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([2]))
    batch_source = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    batch_target = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    batch_dtrans_source = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([d0.shape[2]]))
    batch_dtrans_target = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([d0.shape[2]]))
    found = 0
    while (True):
        sind = np.random.randint(0,nsubjects)
        if len(vol_list[sind]) < 2:
            continue
        tind1 = np.random.randint(0,len(vol_list[sind]))
        tind2 = np.random.randint(0,len(vol_list[sind]))
        while (tind2 == tind1):
            tind2 = np.random.randint(0,len(vol_list[sind]))
        im_src = vol_list[sind][tind1].data
        im_trg = vol_list[sind][tind2].data
        dtrans_src = dtrans_list[sind][tind1]
        dtrans_trg = dtrans_list[sind][tind2]
        if len(im_src.shape) > 2:
            im_src = im_src[...,0]
        if len(im_trg.shape) > 2:
            im_trg = im_trg[...,0]
        batch_source[found,...,0] = im_src / im_src.max()
        batch_target[found,...,0] = im_trg / im_trg.max()
        batch_dtrans_source[found,...] = dtrans_src / dtrans_src.max()
        batch_dtrans_target[found,...] = dtrans_src / dtrans_src.max()
        found += 1
        if found >= batch_size:
            if return_warp == True:
                yield([batch_source, batch_target], [batch_target,warp])
            else:
                if bidir == True:
                    yield([batch_source, batch_target], [batch_target, batch_source])
                else:
                    yield([batch_source, batch_target, batch_dtrans_source, batch_dtrans_target], [batch_target, batch_dtrans_target])
            found = 0


def BatchGeneratorSurf(vol_list, seg_list, dtrans_list, npoints, nlabels, batch_size=4,return_warp=False, bidir=False):
    nsubjects = len(vol_list)
    im0 = vol_list[0][0]
    d0 = dtrans_list[0][0]
    warp = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([2]))

    # targets of mapped points should all be 0
    zpoints = np.zeros(tuple([batch_size, npoints, 1]))
    batch_src = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    batch_trg = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([1]))
    batch_dtrans_src = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([d0.shape[2]]))
    batch_dtrans_trg = np.zeros(tuple([batch_size])+im0.shape[0:2]+tuple([d0.shape[2]]))
    batch_points_src = np.zeros((batch_size,npoints,3))
    batch_points_trg = np.zeros((batch_size,npoints,3))
    found = 0
    while (True):
        sind = np.random.randint(0,nsubjects)
        if len(vol_list[sind]) < 2:
            continue
        tind1 = np.random.randint(0,len(vol_list[sind]))
        tind2 = np.random.randint(0,len(vol_list[sind]))
        while (tind2 == tind1):
            tind2 = np.random.randint(0,len(vol_list[sind]))

        # intensity data
        im_src = vol_list[sind][tind1].data   
        im_trg = vol_list[sind][tind2].data

        # distance transforms
        dtrans_src = dtrans_list[sind][tind1]  
        dtrans_trg = dtrans_list[sind][tind2]

        # segmentations
        seg_src = seg_list[sind][tind1].data
        seg_trg = seg_list[sind][tind2].data

        # find all nonzero indices
        ind_src = np.where(seg_src > 0)
        ind_trg = np.where(seg_trg > 0)

        psrc = np.random.permutation(len(ind_src[0]))
        # handle case where more points than exist in this example src
        while len(psrc) > npoints:  
            np.concatenate(psrc,psrc)
        xp_src = ind_src[0][psrc[0:npoints]]
        yp_src = ind_src[1][psrc[0:npoints]]
        batch_points_src[found,...,0] = xp_src
        batch_points_src[found,...,1] = yp_src
        batch_points_src[found,...,2] = seg_src[xp_src, yp_src]

        ptrg = np.random.permutation(len(ind_trg[0]))
        # handle case where more points than exist in this example trg
        while len(ptrg) > npoints:
            np.concatenate(ptrg,ptrg)
        xp_trg = ind_trg[0][ptrg[0:npoints]]
        yp_trg = ind_trg[1][ptrg[0:npoints]]
        batch_points_trg[found,...,0] = xp_trg
        batch_points_trg[found,...,1] = yp_trg
        batch_points_trg[found,...,2] = seg_trg[xp_trg,yp_trg]
        if len(im_src.shape) > 2:
            im_src = im_src[...,0]
        if len(im_trg.shape) > 2:
            im_trg = im_trg[...,0]
        batch_src[found,...,0] = im_src / im_src.max()
        batch_trg[found,...,0] = im_trg / im_trg.max()
        batch_dtrans_src[found,...] = dtrans_src / dtrans_src.max()
        batch_dtrans_trg[found,...] = dtrans_src / dtrans_src.max()
        found += 1
        if found >= batch_size:
            if return_warp == True:
                inputs = [batch_src, batch_trg]
                outputs = [batch_trg,warp]
            else:
                if bidir == True:
                    inputs = [batch_src, batch_trg]
                    outputs = [batch_trg, batch_src]
                else:  # only this case works at the moment (bidir=False)
                    inputs = [batch_src, batch_trg, batch_dtrans_src, batch_dtrans_trg, batch_points_src, batch_points_trg]
                    outputs = [batch_trg, zpoints,zpoints]
            yield inputs,outputs
            found = 0

train_affine = True
scale=1


batch_size=16

learning_rate = 0.001*.001

transform_type = 'rigid+scale'
transform_type = 'affine'
which_loss = 'tukey'
which_loss = 'ncc'
which_loss = 'mse'
symmetrize = False
affine_blurs = [1]
bidir = False
#model = vxm.networks.VxmAffineDense(target_shape, enc_nf, dec_nf, enc_nf_affine=enc_nf_affine,transform_type=transform_type, affine_bidir=bidir, affine_blurs=affine_blurs)

nb_labels = 2
#model = vxm.networks.VxmAffineSegSemiSupervised(target_shape, enc_nf_affine,transform_type=transform_type, bidir=bidir, blurs=affine_blurs, int_downsize=1, seg_downsize=1, nb_labels=nb_labels)

npoints = 5000
model = vxm.networks.VxmAffineSurfaceSemiSupervised(target_shape, enc_nf_affine,transform_type=transform_type, nb_labels_sample=nb_labels, nb_surface_points=npoints, bidir=bidir, blurs=affine_blurs)

ncc = vxm.losses.NCC([25,25])
sncc = vxm.losses.NCC([25,25])
tukey_c = .25
tukey = vxm.losses.TukeyBiweight(c=tukey_c)
if which_loss == 'ncc':
    aloss=ncc.loss
    sloss=sncc.loss
elif which_loss == 'mse':
    aloss=keras.losses.mse
    sloss=aloss
elif which_loss == 'tukey':
    aloss = tukey.loss
    sloss=aloss

ldir = 'logs/cxr.%s.train_affine.%s' % (which_loss, str(train_affine))
tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=ldir)

reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='loss', factor=0.5, patience=10, min_lr=1e-8, verbose=True, min_delta=0.0001)
affine_wt_fname = 'cxr.vxm.mpr.%s.%s.sym.%s.train_affine.%s.h5' % (transform_type, which_loss, str(symmetrize),str(train_affine))
callbacks = [ModelSaver(model, 100, affine_wt_fname,cp_iters=10), reduce_lr, tensorboard_callback]
steps_per_epoch = int(128)

#affine_model = model.references.affine_model
affine_model = model

nlabels = 2 # not included 0 label
bg = BatchGeneratorSurf(vol_list, seg_list, dtrans_list, npoints, nlabels, batch_size=batch_size,return_warp=False, bidir=bidir)

losses = [aloss,aloss,aloss]
affine_model.compile(optimizer=keras.optimizers.Adam(lr=.1*learning_rate), loss=losses,loss_weights=[.00001, 1,1])
initial_epoch = 0
epochs = initial_epoch+50
fithistr = affine_model.fit_generator(bg, steps_per_epoch = steps_per_epoch, epochs = epochs, verbose=1, callbacks=callbacks, class_weight=None)
if 1:
    np.random.seed(1)
    bg2 = BatchGeneratorSurf(vol_list, seg_list, dtrans_list, npoints, nlabels, batch_size=batch_size,return_warp=False, bidir=bidir)
    inb,outb = next(bg2)
    p = affine_model.predict(inb)
    fv = sf.vis.Freeview()
    fv.add_image(np.transpose(inb[0][...,0], (1,2,0)), name='src')
    fv.add_image(np.transpose(inb[1][...,0], (1,2,0)), name='trg')
    if isinstance(p, list):
        p = p[0]
    fv.add_image(np.transpose(p[...,0], (1,2,0)), name='p')
    fv.show()
    affine_pred_model = keras.models.Model(affine_model.inputs, affine_model.references.affines[0])
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

nl_wt_fname = 'cxr.vxm.mpr.%s.%s.sym.%s.train_affine.%s.nl.h5' % (transform_type, which_loss, str(symmetrize),str(train_affine))
which_loss = 'mse'
aloss=keras.losses.mse
losses = [aloss, vxm.losses.Grad().loss]
smoothness_weight=.001*100
model.compile(optimizer=keras.optimizers.Adam(lr=.1*learning_rate), loss=losses, loss_weights=[1,smoothness_weight])

bg_nl = BatchGenerator(vol_list, dtrans_list, batch_size=batch_size,return_warp=True, bidir=False)

nl_epochs = 50
initial_epoch = 0
epochs = initial_epoch+nl_epochs
fithistr = model.fit_generator(bg_nl, steps_per_epoch = steps_per_epoch, epochs = initial_epoch+epochs, verbose=1, callbacks=callbacks, class_weight=None, initial_epoch=initial_epoch)
affine_model.save(affine_wt_fname)
model.save(nl_wt_fname)
if 1:
    np.random.seed(1)
    bg_nl2 = BatchGeneratorSurf(vol_list, seg_list, dtrans_list, npoints, nlabels, batch_size=batch_size,return_warp=False, bidir=bidir)

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
        name = 'pred.%s.%s.sym.%s.affine.%s' % (transform_type,str(which_loss)+str(tukey_c),str(symmetrize),str(train_affine))
    else:
        name = 'pred.%s.%s.sym.%s.affine.%s' % (transform_type,str(which_loss),str(symmetrize),str(train_affine))
    fv.add_image(np.transpose(affine_p, (1,2,0,3)), name='affine'+name)
    fv.add_image(np.transpose(nl_p, (1,2,0,3)), name='NL'+name)
    fv.show()
    affine_pred_model = affine_model.get_predictor_model()
    affine = affine_pred_model.predict(inb)
    maffine = np.reshape(np.append(affine[0,:],np.array([0,0,1])),(3,3))
    pprint(maffine)
    w,v = np.linalg.eig(maffine)
    print(w,v)
