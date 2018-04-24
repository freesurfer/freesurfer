import pdb as gdb
import tensorflow as tf
from keras import backend as K
import nibabel as nib
from nibabel import processing as nip
import unet_model
from unet_model import class_net
import numpy as np
import matplotlib.pyplot as plt
import keras
import os

sdir = os.getenv('SUBJECTS_DIR')
slist = os.getenv('SURF_AUTO_LABEL').split(' ')
nsubjects = np.array(slist).shape[0]
#nsubjects = 10

for sno in range(nsubjects):
    subject = '%s.gentle' % slist[sno]
    print 'processing subject %d of %d: %s' % (sno, nsubjects, subject.strip())
    pdir = '%s/%s/patches' % (sdir, subject)
    rh_patches = nib.load('%s/rh.patches.mgz' % pdir)
    rh_labels = nib.load('%s/rh.labels.mgz' % pdir)
    lh_patches = nib.load('%s/lh.patches.mgz' % pdir)
    lh_labels = nib.load('%s/lh.labels.mgz' % pdir)
    print 'subject has %d patches' % rh_patches.get_data().shape[3]
    if (sno == 0):
        patches = lh_patches.get_data() 
        labels = lh_labels.get_data()[:,1,:]
    else:
        patches = np.concatenate((patches, lh_patches.get_data()), axis=3)
        labels = np.concatenate((labels, lh_labels.get_data()[:,1,:]))
    patches = np.concatenate((patches, rh_patches.get_data()), axis=3)
    labels = np.concatenate((labels, rh_labels.get_data()[:,1,:]))
    print 'appending patches, total now %d' % patches.shape[3]

learning_rate = 0.00001 ;
nfilters = 30;
ndim = patches.ndim-1
patch_shape = (patches.shape[0], patches.shape[1], patches.shape[2])
net = class_net(patch_shape, ndim, nfilters, 2, learning_rate)

reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='loss', factor=0.5 ,patience=5, min_lr=0.001*learning_rate)
x_train = np.transpose(patches, (3, 0, 1, 2));
x_train = np.reshape(x_train, (x_train.shape[0], x_train.shape[1], x_train.shape[2], x_train.shape[3], 1))

fit = net.fit(x=x_train, y=labels, batch_size=20, epochs=200, verbose=1, callbacks=[reduce_lr],validation_split=.1)
net.save_weights('fcd.net')
