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
import imageio, pydicom, gdcm, load_serial_cxr

bdir = '/autofs/cluster/lcnextdata1/CCDS_CXR/CXR-Serial/def_20200413'


il, sl, sn = load_serial_cxr.load_serial_cxr(bdir)

for sno, ilist in enumerate(il):
    if len(ilist)>=2:
        date_list = []
        time_list = []
        for ino, im in enumerate(ilist):
            date_list.append(int(im.StudyDate))
            if hasattr(im, 'SeriesTime'):
                time_list.append(int(im.SeriesTime.split('.')[0]))
            else:
                time_list.append(int(im.StudyTime.split('.')[0]))

        ind = np.array(date_list).argsort()
        date_list2 = []
        time_list2 = []
        ilist2 = []
        for i in ind.tolist():
            date_list2.append(date_list[i])
            time_list2.append(time_list[i])
            ilist2.append(ilist[i])

        for ino, im in enumerate(ilist2):
            tokens = sl[sno][ind[ino]].split('/')
            fname = '/'.join(tokens[0:-2]) + '/time%2.2d.mgz' % ino
            sf.Volume(im.pixel_array.astype(np.float32)).save(fname)
