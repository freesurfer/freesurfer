"""
Test models for MICCAI 2018 submission of VoxelMorph.
"""

# py imports
import os
import sys
import glob

# third party
import tensorflow as tf
import scipy.io as sio
import numpy as np
import keras
from keras.backend.tensorflow_backend import set_session
from scipy.interpolate import interpn

# project
sys.path.append('../ext/medipy-lib')
import medipy
import network2d
# import util
from medipy.metrics import dice
import datagenerators



def test(gpu_id, model_dir, iter_num, data_dir, file_name,
         compute_type = 'GPU',  # GPU or CPU
         vol_size=(528,256),
         nf_enc=[16,32,32,32],
         nf_dec=[32,32,32,32,16,3],
         save_file=None):
    """
    test via segmetnation propagation
    works by iterating over some iamge files, registering them to atlas,
    propagating the warps, then computing Dice with atlas segmentations
    """

    # atlas files
    atlas = sio.loadmat('../atlasdata/lh.pad.mat')
    atlas_vol = atlas['tmp'][np.newaxis,...,np.newaxis]

    # GPU handling
    gpu = '/gpu:' + str(gpu_id)
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.allow_soft_placement = True
    set_session(tf.Session(config=config))

    # load weights of model
    with tf.device(gpu):
        # if testing miccai run, should be xy indexing.
        net = network2d.miccai2018_net(vol_size, nf_enc,nf_dec, use_miccai_int=True, indexing='xy')
        net.load_weights(os.path.join(model_dir, str(iter_num) + '.h5'))

        # compose diffeomorphic flow output model
        diff_net = keras.models.Model(net.inputs, net.get_layer('diffflow').output)

        # NN transfer model
        nn_trf_model = network2d.nn_trf(vol_size)

    # if CPU, prepare grid
    if compute_type == 'CPU':
        grid, xx, yy, zz = util.volshape2grid_3d(vol_size, nargout=4)



    D = datagenerators.load_surf_by_name(data_dir + file_name)

    #X_vol = D[:,:,:,1,:]
    X_vol = D


    with tf.device(gpu):
        pred = diff_net.predict([ X_vol, atlas_vol])

    print(pred.shape)

    # Warp segments with flow
    if compute_type == 'CPU':
        flow = pred[0, :, :, :]
        warp_seg = util.warp_seg(X_vol, flow, grid=grid, xx=xx, yy=yy)

    else:  # GPU
        warp_seg = nn_trf_model.predict([X_vol, pred])[0,...,0]

    print(warp_seg.shape)

    flow = pred[0, :, :, :]
    sio.savemat('flow.mat',{'flow':flow})

    sio.savemat('y.mat',{"tmp":warp_seg})


if __name__ == "__main__":
    """
    assuming the model is model_dir/iter_num.h5
    python test_jc.py gpu_id model_dir data_dir iter_num
    """
    test(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
