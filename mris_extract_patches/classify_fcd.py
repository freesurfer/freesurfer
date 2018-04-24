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

sd = os.getenv('SUBJECTS_DIR')
slist = os.getenv('SURF_AUTO_LABEL').split(' ')
nsubjects = np.array(slist).shape[0]
learning_rate = 0.00001 ;
nfilters = 30;
batch_size = 10

for sno in range (7,nsubjects):
    with tf.Session() as sess:
        subject = '%s.gentle' % slist[sno]
        print 'processing subject %d of %d: %s' % (sno, nsubjects, subject)
        pdir = '%s/%s/cortex_patches' % (sd, subject)
        sdir = '%s/%s/surf' % (sd, subject)

        hemi = 'lh'
        lh_patches = nib.load('%s/%s.patches.mgz' % (pdir,hemi))
        lh_labels = nib.load('%s/%s.labels.mgz' % (pdir,hemi))
        print 'subject has %d patches' % lh_patches.get_data().shape[3]

        patches = lh_patches.get_data() 
        vertices = lh_labels.get_data()[:,0,:]
        labels = lh_labels.get_data()[:,1,:]
        print 'appending patches, total now %d' % patches.shape[3]

        ndim = patches.ndim-1
        patch_shape = (patches.shape[0], patches.shape[1], patches.shape[2])
        net = class_net(patch_shape, ndim, nfilters, 2, learning_rate)
        net.load_weights('fcd.net')

        x_input = np.transpose(patches, (3, 0, 1, 2));
        x_input = np.reshape(x_input, (x_input.shape[0], x_input.shape[1], x_input.shape[2], x_input.shape[3], 1))
        nvertices = x_input.shape[0]

        surface = nib.freesurfer.io.read_geometry('%s/%s.white' % (sdir,hemi))
        surf_vertices = surface[0].shape[0]
        
        pvals = np.zeros((surf_vertices,1))
        sulc = nib.freesurfer.io.read_morph_data('%s/%s.sulc' % (sdir,hemi))
        for vno in range(0, nvertices, batch_size):
            if ((vno % 5000) == 0):
                print 'processing vno %d of %d: %2.2f%%' % (vno, nvertices, 100.0*vno/nvertices)
            if sulc[vertices[vno]] < 0:
                continue
            pvals_batch = net.predict(x_input[vno:vno+batch_size,:,:,:,:], batch_size, verbose=0)
            batch_vertices = vertices[vno:vno+batch_size]
            pvals[batch_vertices] = np.reshape(pvals_batch, (pvals_batch.shape[0],1,1))

        surface_overlay = nib.Nifti1Image(np.reshape(pvals,(surf_vertices,1,1)),np.diag([1, 1, 1, 1]));
        nib.save(surface_overlay, '%s/%s.fcd_overlay.mgz' % (pdir,hemi))
        
        hemi = 'rh'
        rh_patches = nib.load('%s/%s.patches.mgz' % (pdir,hemi))
        rh_labels = nib.load('%s/%s.labels.mgz' % (pdir,hemi))
        patches = rh_patches.get_data() 
        vertices = rh_labels.get_data()[:,0,:]
        labels = rh_labels.get_data()[:,1,:]
        x_input = np.transpose(patches, (3, 0, 1, 2));
        x_input = np.reshape(x_input, (x_input.shape[0], x_input.shape[1], x_input.shape[2], x_input.shape[3], 1))
        
        surface = nib.freesurfer.io.read_geometry('%s/%s.white' % (sdir,hemi))
        surf_vertices = surface[0].shape[0]
    
        pvals = np.zeros((surf_vertices,1))
        sulc = nib.freesurfer.io.read_morph_data('%s/%s.sulc' % (sdir,hemi))
        for vno in range(0, nvertices, batch_size):
            if ((vno % 5000) == 0):
                print 'processing vno %d of %d: %2.2f%%' % (vno, nvertices, 100.0*vno/nvertices)
            if sulc[vertices[vno]] < 0:
                continue
            pvals_batch = net.predict(x_input[vno:vno+batch_size,:,:,:,:], batch_size, verbose=0)
            batch_vertices = vertices[vno:vno+batch_size]
            pvals[batch_vertices] = np.reshape(pvals_batch, (pvals_batch.shape[0],1,1))
    
        surface_overlay = nib.Nifti1Image(np.reshape(pvals,(surf_vertices,1,1)),np.diag([1, 1, 1, 1]));
        nib.save(surface_overlay, '%s/%s.fcd_overlay.mgz' % (pdir,hemi))

