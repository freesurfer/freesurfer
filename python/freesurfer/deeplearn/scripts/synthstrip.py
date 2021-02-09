
"""
Example script for testing quality of trained vxm models. This script iterates over a list of images
and corresponding segmentations, registers them to an atlas, propagates segmentations to the atlas,
and computes the dice overlap. Example usage is:

    test.py \
    --model models/model.h5 \
    --invol data/withskull.mgz \

Where each atlas and scan npz file is assumed to contain the array variables 'vol' and 'seg'. This
script will most likely need to be customized to fit your data.

If you use this code, please cite the following, and read function docs for further info/citations
    VoxelMorph: A Learning Framework for Deformable Medical Image Registration 
    G. Balakrishnan, A. Zhao, M. R. Sabuncu, J. Guttag, A.V. Dalca. 
    IEEE TMI: Transactions on Medical Imaging. 38(8). pp 1788-1800. 2019. 

Copyright 2020 Adrian V. Dalca

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in
compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing permissions and limitations under the
License.
"""

import os
import argparse
import numpy as np
import voxelmorph as vxm
import tensorflow as tf
import copy

from freesurfer import deeplearn as fsd
import freesurfer as fs

import neurite as ne
import neurite_sandbox as nes
import voxelmorph as vxm
import read_long
from netparms import *
from tensorflow.keras import layers as KL


# parse commandline args
parser = argparse.ArgumentParser()
parser.add_argument('--model', required=True, help='keras model file')
parser.add_argument('--wts', help='wt filename')
parser.add_argument('--invol', required=True, help='input MRI volume')
parser.add_argument('--outvol', required=True, help='output MRI volume')
parser.add_argument('--pred', help='write prediction volume')
parser.add_argument('--norm', help='write normalized volume')
parser.add_argument('--uthresh', type=float, help='specify threshold to erase above')
parser.add_argument('--border', help='number of border voxels to set threshold at', default=4,type=int)
parser.add_argument('--gpu', help='GPU number - if not supplied, CPU is used')
parser.add_argument('--multichannel', action='store_true',
                    help='specify that data has multiple channels')
args = parser.parse_args()

# read input volume and normalize input intensities
mri_in = fs.Volume.read(args.invol)
if args.uthresh:
    mri_in.data[mri_in.data>args.uthresh] = 0

in_data = (mri_in.data - mri_in.data.min()) 
in_data = np.clip(in_data / np.percentile(in_data, 97), 0,1)
mri_in.data = in_data
if 0:
    mri_conf = mri_in.reslice((1,1,1))
    left_pad = (np.array((256,256,256))-np.array(mri_conf.shape))//2
    right_pad = (np.array((256,256,256))-(np.array(mri_conf.shape)+left_pad))
    padding = ((left_pad[0], right_pad[0]), (left_pad[1], right_pad[1]), (left_pad[2], right_pad[2]))
    conf_data = np.pad(mri_conf.data, padding)
    #mri_conf = mri_conf1.conform_to_shape((256,256,256))
else:
    conf_data = mri_in.data
    mri_conf = mri_in



if args.norm:
    mri_in.write('norm.mgz')

# device handling, model reading and prediction
if args.gpu:
    device, ngpus = vxm.tf.utils.setup_device(args.gpu)
else:
    device = '/cpu:0'

print(f'using device {device}')
with tf.device(device):
    print(f'loading model from {args.model}')
    synthmodel = ne.models.SynthStrip.load(args.model)
    model = synthmodel.get_strip_model()
    if args.wts:
        print(f'loading weights from {args.wts}')

    synthmodel.load_weights(args.wts)
    pred = model.predict(conf_data[np.newaxis,...,np.newaxis]).squeeze()

if 0:
    unconf_data = pred[left_pad[0]:-right_pad[0], left_pad[1]:-right_pad[1], left_pad[2]:-right_pad[2]]
    mri_unconf1 = copy.copy(mri_conf)
    mri_unconf1.data = unconf_data
    mri_unconf = mri_unconf1.reslice(mri_in.voxsize).conform_to_shape(mri_in.shape)
else:
    mri_mask = copy.copy(mri_conf)
    mri_mask.data = pred


# create and write output volume
mri_out = copy.deepcopy(mri_in)
mri_out.data = mri_in.data * (mri_mask.data.squeeze() < args.border).astype(np.float)
mri_out.write(args.outvol)

if args.pred:
    mri_mask.write(args.pred)
