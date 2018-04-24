import os
import glob
import numpy as np
import tables
import nibabel as nib
import tempfile
from random import shuffle

import sys
from os.path import join as opj
from deeplearn_utils import DeepImageSynth

import socket
import subprocess


def robust_normalize(in_img_data):
    in_img_flat = np.ravel(in_img_data, 'C')
    in_img_fg = in_img_flat[in_img_flat > 0].reshape(-1, 1)
    p99 = np.percentile(in_img_fg, q=99)

    # set p99 to 255
    scaling = 255.0 / p99
    print(scaling)
    out_img_data = in_img_data * scaling
    return out_img_data


if __name__ == "__main__":
    import os

    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  # see issue #152
    os.environ["CUDA_VISIBLE_DEVICES"] = "2"

    training_dirs = ['/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20171129_hu006_mpr_t2sp_fastmov',
                   '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20171106_trac_hu_mprage002',
                   '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20171113_trac_human_mprage003_slowmov',
                   '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20171114_trac_hu_mprage004',
                   '/autofs/space/vault_020/users/rfrost/moco/tracoline/rawdata/20171130_mpr_epi']

    still_filenames = []
    mov_filenames = []
    moco_filenames = []

    for tdir in training_dirs:
        still_filenames = still_filenames + sorted(glob.glob(os.path.join(tdir, '*mprage*still*.nii.gz')))
        mov_filenames = mov_filenames + sorted(glob.glob(os.path.join(tdir, '*mprage*mov_mocoOFF*.nii.gz')))
        moco_filenames = moco_filenames + sorted(glob.glob(os.path.join(tdir, '*mprage*mov*mocoON*.nii.gz')))


    all_filenames = still_filenames + mov_filenames
    motion_labels = np.vstack(( np.zeros((len(still_filenames),1)), np.ones((len(mov_filenames),1))  ))
    feature_shape = (256, 256)
    training_idxs = np.array([0, 1, 2, 3, 6, 7, 16, 17, 18, 19,20,21])
    validation_idxs = np.array([6,7,20,21])
    testing_idxs = np.array([8,9,10,11,12,13,14,15,22,23,24,25])



    training_filenames = list()
    validation_filenames = list()
    testing_filenames = list()

    for idx in training_idxs:
        training_filenames.append(all_filenames[idx])

    for idx in validation_idxs:
        validation_filenames.append(all_filenames[idx])

    for idx in testing_idxs:
        testing_filenames.append(all_filenames[idx])


    d  = 6
    f = 32
    hostname = socket.gethostname()
    if  hostname == 'bhim.nmr.mgh.harvard.edu':
        tmp_folder = '/local_mount/space/bhim/1/users/aj660/tmp'
    elif hostname == 'nike.nmr.mgh.harvard.edu':
        tmp_folder = '/local_mount/space/nike/1/users/amod/tmp'
    elif hostname == 'trabant.nmr.mgh.harvard.edu':
        tmp_folder = '/tmp/aj660'
        subprocess.call(['mkdir', '-p', tmp_folder])


    curr_unet = DeepImageSynth.DeepImageSynth(net='class_net',
                                              dim=2,
                                              unet_num_filters=f, unet_depth=d,
                                              unet_downsampling_factor=1,
                                              feature_shape=feature_shape,
                                              storage_loc="memory",
                                              temp_folder=tmp_folder,
                                              n_labels=2,
                                              loss='binary_crossentropy',
                                              use_patches=False,
                                              wmp_standardize=True,
                                              fcn=False,
                                              )
    curr_unet.load_training_slices_and_labels(training_filenames,
                                              label_list=motion_labels[training_idxs],
                                              is_src_label_img=False)
    curr_unet.load_validation_slices_and_labels(validation_filenames,
                                                label_list=motion_labels[validation_idxs],
                                                is_src_label_img=False)

    output_dir = "/autofs/space/mreuter/users/amod/deep_learn/results/motion_detect_"+ \
                 str(d) +"_filts_" + str(f) +"_patch"+str(feature_shape[0])+"x"+str(feature_shape[1])
    subprocess.call(['mkdir', '-p', output_dir])
    print(curr_unet.model.summary())
    curr_unet.train_network(output_prefix=opj(output_dir, 'motion_detect_'), epochs=100,
                            initial_epoch=1,batch_size=32, steps_per_epoch=10000,
                            save_per_epoch=True, save_weights=True)


    model_dir = '/autofs/space/mreuter/users/amod/deep_learn/results/motion_detect_6_filts_32_patch256x256'


    # net_file = opj(model_dir, 'motion_detect__model_epoch3.h5')
    # loss = 'binary_cross_entropy'
    # curr_unet = DeepImageSynth.DeepImageSynth.from_file_old_model(net_file,loss, n_labels=2,
    #                                                     storage_loc='disk', temp_folder=model_dir)

    predicted_labels = list()
    true_labels = list()
    predicted_probs = list()
    for test_idx, in_file in enumerate(testing_filenames):
        print(in_file)
        curr_pred_probs, curr_pred_labels = curr_unet.predict_labels(in_file)
        predicted_labels.append(curr_pred_labels)
        predicted_probs.append(curr_pred_probs)

        curr_true_label = motion_labels[testing_idxs[test_idx]]
        curr_true_labels = np.zeros(curr_pred_labels.shape)
        curr_true_labels = curr_true_labels * curr_true_label
        true_labels.append(curr_true_labels)


    predicted_labels = np.concatenate(predicted_labels)
    true_labels = np.concatenate(true_labels)

    true_positives = len(np.where((predicted_labels == 1) & (true_labels == 1))[0])
    true_negatives = len(np.where((predicted_labels == 0) & (true_labels == 0))[0])
    false_positives = len(np.where((predicted_labels == 1) & (true_labels == 0))[0])
    false_negatives = len(np.where((predicted_labels == 0) & (true_labels == 1))[0])

    confusion_matrix = np.zeros((2,2))
    confusion_matrix[0, 0] = true_negatives
    confusion_matrix[1, 0] = false_positives
    confusion_matrix[0, 1] = false_negatives
    confusion_matrix[1, 1] = true_positives

    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)

