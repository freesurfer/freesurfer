"""
train atlas-based alignment with MICCAI2018 version of VoxelMorph,
Input training spherical images were after padding along longitude
"""

# python imports
import os
import glob
import sys
import random
from argparse import ArgumentParser

# third-party imports
import tensorflow as tf
import numpy as np
from keras.backend.tensorflow_backend import set_session
from keras.optimizers import Adam
from keras.models import load_model, Model
import scipy.io as sio

# project imports
import datagenerators
import network2d
import losspad


def train(data_dir, model_dir, gpu_id, lr, n_iterations, alpha, model_save_iter, gamma=10000, batch_size=1):
    """
    model training function
    :param model_dir: model folder to save to
    :param gpu_id: integer specifying the gpu to use
    :param lr: learning rate
    :param n_iterations: number of training iterations
    :param alpha: the alpha, the scalar in front of the smoothing laplacian
    :param model_save_iter: frequency with which to save models
    :param batch_size: Optional, default of 1. can be larger, depends on GPU memory and volume size
    """

    # load Buckner atlas from provided files. This atlas is 528*256.
    atlas1 = sio.loadmat('../atlasdata/lh.pad.mat')
    atlas_vol1 = atlas1['tmp'][np.newaxis,...,np.newaxis]


    sigma1 = sio.loadmat('../atlasdata/weight.pad.mat')
    image_sigma1 = sigma1['tmp'][np.newaxis,...,np.newaxis]


    vol_size = atlas_vol1.shape[1:-1]

    # prepare model folder
    if not os.path.isdir(model_dir):
        os.mkdir(model_dir)
    #print(model_dir)

    # gpu handling
    gpu = '/gpu:' + str(gpu_id)
    os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu_id)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.allow_soft_placement = True
    set_session(tf.Session(config=config))

    train_surf_names = glob.glob(data_dir + '/**/*.sphere')
    random.shuffle(train_surf_names)  # shuffle volume list

    # Diffeomorphic network architecture used in MICCAI 2018 paper
    nf_enc = [16,32,32,32]
    nf_dec = [32,32,32,32,16,3]

    # prepare the model
    # in the CVPR layout, the model takes in [image_1, image_2] and outputs [warped_image_1, velocity_stats]
    # in the experiments, we use image_2 as atlas
    with tf.device(gpu):
        # miccai 2018 used xy indexing.
        model = network2d.miccai2018_net(vol_size,nf_enc,nf_dec, use_miccai_int=True, indexing='xy')

        # compile
        model_losses = [losspad.kl_l2loss(image_sigma1), losspad.kl_loss(alpha), losspad.bound_loss(gamma)]

	#model_losses = loss2d.kl_l2loss(image_sigma)
        model.compile(optimizer=Adam(lr=lr), loss=model_losses)

        #model.load_weights('../hcpmodel/30000.h5')

        # save first iteration
        model.save(os.path.join(model_dir,  str(0) + '.h5'))

    #print(model.summary())

    train_example_gen = datagenerators.sphere_gen(train_surf_names)
    zeros = np.zeros((1, *vol_size, 2))

    # train. Note: we use train_on_batch and design out own print function as this has enabled
    # faster development and debugging, but one could also use fit_generator and Keras callbacks.
    for step in range(1, n_iterations):

        # get_data
        D = next(train_example_gen)[0]

        #pred = model.predict([X,atlas_vol])

        # train
        with tf.device(gpu):

            train_loss = model.train_on_batch([D,atlas_vol1], [atlas_vol1, zeros, zeros])

        if not isinstance(train_loss,list):
            train_loss = [train_loss]

        # print
        print_loss(step, 0, train_loss)

        # save model
        with tf.device(gpu):
            if (step % model_save_iter == 0) or step < 10:
                model.save(os.path.join(model_dir,  str(step) + '.h5'))


def print_loss(step, training, train_loss):
    """
    Prints training progress to std. out
    :param step: iteration number
    :param training: a 0/1 indicating training/testing
    :param train_loss: model loss at current iteration
    """
    s = str(step) + "," + str(training)

    if isinstance(train_loss, list) or isinstance(train_loss, np.ndarray):
        for i in range(len(train_loss)):
            s += "," + str(train_loss[i])
    else:
        s += "," + str(train_loss)

    print(s)
    sys.stdout.flush()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--model_dir", type=str,
                        dest="model_dir", default='../model/',
                        help="models folder")
    parser.add_argument("--gamma", type=float, default=10000,
                        dest="gamma", help="boundary term weight")
    parser.add_argument("--gpu", type=int, default=0,
                        dest="gpu_id", help="gpu id number")
    parser.add_argument("--lr", type=float,
                        dest="lr", default=1e-4, help="learning rate")
    parser.add_argument("--iters", type=int,
                        dest="n_iterations", default=150001,
                        help="number of iterations")
    parser.add_argument("--alpha", type=float,
                        dest="alpha", default=30000000,
                        help="alpha regularization parameter")
    parser.add_argument("--data_dir", type=str,
                        dest="data_dir", default='traindata/',
                        help="training data folder")
    parser.add_argument("--checkpoint_iter", type=int,
                        dest="model_save_iter", default=10000,
                        help="frequency of model saves")

    args = parser.parse_args()
    train(**vars(args))
