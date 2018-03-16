import os

import glob
import numpy as np
import tables
import nibabel as nib
import tempfile
from random import shuffle

import sys
from os.path import join as opj

import subprocess
from keras import backend
from keras.models import load_model
from keras.optimizers import serialize
from keras.callbacks import Callback
from deeplearn_utils import DeepImageSynth
from nipype.interfaces.base import split_filename
import pandas as pd
import simplejson
import time
from nipype.interfaces.base import split_filename


def fetch_training_data_files(subjects_dir, img_input_type, training_subject_idxs):
    ''' assumes a freesurfer directory structure
    # Arguments
    :param fs_dir: directory with all the scanner freesurfer results stored
    :param src_scanner: scanner directory name from which we want to extract training images
    :param: src_img_input_type: freesurfer output that we want to extract for e.g. orig/001.mgz or aseg.mgz etc.
    :param trg_scanner
    :param trg_img_input_types tuple('orig/001', 'aparc+aseg') etc.
    '''
    training_data_files = list()
    # input_subj_dir_list = list()
    # for i in training_subject_idxs:
    #     input_subj_dir_list.append(src_subj_dir_list[i])
    for training_subject_id in training_subject_idxs:
        training_data_files.append(os.path.join(subjects_dir,training_subject_id,'mri', img_input_type + ".mgz"))
    return training_data_files



if __name__ == "__main__":

    # get all the aseg labels

    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  # see issue #152
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    aseg_labels =np.loadtxt('aseg_labels.txt')

    subjects_dir = '/autofs/space/mreuter/users/amod/deep_learn/data/aseg_atlas'

    # get enough young, mid, old and AD brains. Half? load the demographics file
    demographics_file = opj(subjects_dir, 'demographics')
    demo_df = pd.read_csv(demographics_file, index_col=0, delim_whitespace=True)
    subject_id_list = demo_df.index.get_values()

    young_ids = demo_df.index[demo_df['group']=='YOUNG']
    mid_ids = demo_df.index[demo_df['group']=='MID']
    old_ids = demo_df.index[demo_df['group']=='OLD']
    ad1_ids = demo_df.index[demo_df['group']=='AD_1']
    adp5_ids = demo_df.index[demo_df['group'] == 'AD_p5']

    np.random.seed(1729)
    y_idxs = sorted(np.random.choice(len(young_ids), size=(5,), replace=False))
    m_idxs = sorted(np.random.choice(len(mid_ids), size=(5,), replace=False))
    o_idxs = sorted(np.random.choice(len(old_ids), size=(4,), replace=False))
    a1_idxs = sorted(np.random.choice(len(ad1_ids), size=(2,), replace=False))
    adp5_idxs = sorted(np.random.choice(len(adp5_ids), size=(3,), replace=False))

    young_training_ids = young_ids[y_idxs]
    mid_training_ids = mid_ids[m_idxs]
    old_training_ids = old_ids[o_idxs]
    ad1_training_ids = ad1_ids[a1_idxs]
    adp5_training_ids = adp5_ids[adp5_idxs]

    training_id_list = list(young_training_ids[:-1]) + list(mid_training_ids[:-1]) + list(old_training_ids[:-1])\
                       + list(ad1_training_ids) + list(adp5_training_ids)

    validation_id_list = [young_training_ids[-1]] +[mid_training_ids[-1]] + [(old_training_ids[-1])]

    src_img_input_type = 'nu_masked'
    src_filenames = fetch_training_data_files(subjects_dir, src_img_input_type, training_id_list)

    trg_img_input_type = 'seg_edited'
    trg_filenames = fetch_training_data_files(subjects_dir, trg_img_input_type, training_id_list)


    print('src_filenames')
    print(src_filenames)
    print('trg_filenames')
    print(trg_filenames)


    val_src_filenames = fetch_training_data_files(subjects_dir, src_img_input_type, validation_id_list)
    val_trg_filenames = fetch_training_data_files(subjects_dir, trg_img_input_type, validation_id_list)

     # aparcaseg_labels = np.loadtxt('aparc+aseg_labels.txt')

    feature_shape = (32, 32, 32)
    f = 32
    d = 4
    curr_unet = DeepImageSynth.DeepImageSynth(unet_num_filters=f, unet_depth=d,
                                              unet_downsampling_factor=1,
                                              feature_shape=feature_shape,
                                              storage_loc="memory",
                                              temp_folder="/autofs/space/bhim_001/users/aj660/tmp",
                                              loss='dice_coef_loss2',
                                              initial_learning_rate=0.001,
                                              n_labels=len(aseg_labels),
                                              labels=aseg_labels,
                                              num_gpus=1)

    output_dir = "/autofs/space/mreuter/users/amod/cnn-aseg-demo/results/demo_aseg_atlas_unet_depth_" + \
                 str(d) + "_filts_" + str(f) + "_patch" + str(feature_shape[0]) + "x" + str(feature_shape[1]) + \
                 "x" + str(feature_shape[2])
    subprocess.call(['mkdir', '-p', output_dir])


    step_size = (4, 4, 4)

    curr_unet.load_training_images(src_filenames, trg_filenames, is_src_label_img=False,
                                   is_trg_label_img=True, step_size=step_size)

    curr_unet.load_validation_images(val_src_filenames, val_trg_filenames,
                                     is_src_label_img=False,
                                     is_trg_label_img=True,
                                     step_size=step_size)



    # save training validation and testing info
    curr_time = time.strftime("%Y-%m-%d-%H-%M")
    f = open(opj(output_dir, 'training_info_'+ curr_time + '.txt'), 'w')
    simplejson.dump(training_id_list, f)
    f.close()

    f = open(opj(output_dir, 'validation_info_'+ curr_time + '.txt'), 'w')
    simplejson.dump(validation_id_list, f)
    f.close()


    print(curr_unet.model.summary())
    curr_unet.train_network(output_prefix=opj(output_dir, 'unet_aseg_network_'+curr_time), epochs=100,
                            initial_epoch=1,batch_size=32, steps_per_epoch=10000, validation_steps=1000,
                            save_per_epoch=True, save_weights=True)

    #
    # model_dir = '/autofs/space/mreuter/users/amod/deep_learn/results/' \
    #             'seg5_unet_depth_4_filts_32_patch32x32x32/' \

    #
    # net_file = opj(model_dir, 'TRIOmechoandTRIOmprage_1_to_TRIOmprage_1_model_epoch4.h5')
    # loss = 'dice_coef_loss2'
    # curr_unet = DeepImageSynth.DeepImageSynth.from_file(net_file,loss, n_labels=6,
    #                                                     storage_loc='disk', temp_folder=model_dir)
    testing_id_list = list(set(subject_id_list).difference(set(training_id_list)))
    test_src_filenames = fetch_training_data_files(subjects_dir, src_img_input_type, testing_id_list)


    for test_file in test_src_filenames:
        subj_id = split_filename(test_file)[0].split("/")[-2]  # subject id location
        print("Processing " + subj_id)
        out_dir = opj(output_dir, subj_id)
        subprocess.call(['mkdir', '-p', out_dir])
        out_membership_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_soft.mgz")
        out_hard_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_hard.mgz")

        curr_unet.predict_segmentation(test_file, out_membership_file, out_hard_file, step_size=[28, 28, 28])

    for test_file in val_src_filenames:
        subj_id = split_filename(test_file)[0].split("/")[-2]  # subject id location
        print("Processing " + subj_id)
        out_dir = opj(output_dir, subj_id)
        subprocess.call(['mkdir', '-p', out_dir])
        out_membership_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_soft.mgz")
        out_hard_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_hard.mgz")

        curr_unet.predict_segmentation(test_file, out_membership_file, out_hard_file, step_size=[28, 28, 28])

    for test_file in src_filenames:
        subj_id = split_filename(test_file)[0].split("/")[-2]  # subject id location
        print("Processing " + subj_id)
        out_dir = opj(output_dir, subj_id)
        subprocess.call(['mkdir', '-p', out_dir])
        out_membership_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_soft.mgz")
        out_hard_file = opj(out_dir, split_filename(test_file)[1] + "_aseg_unet_depth_" + str(d)
                       + "_filts_" + str(f) + "_hard.mgz")

        curr_unet.predict_segmentation(test_file, out_membership_file, out_hard_file, step_size=[28, 28, 28])
#
    # for test_file in test_src_filenames:
    #     subj_id = split_filename(test_file)[0].split("/")[-3]  # subject id location
    #     print("Processing " + subj_id)
    #     out_dir = opj(output_dir, src_scanner, subj_id)
    #     subprocess.call(['mkdir', '-p', out_dir])
    #     out_file = opj(out_dir, split_filename(test_file)[1] + "_unet_depth_" + str(d)
    #                    + "_filts_" + str(f) + "_trainid_0_2_8_9.mgz")
    #     curr_unet.synthesize_image(test_file, out_file, step_size=[4, 4, 4])
    #
    # trg_img_input_type = 'orig/001'
    # trg_scanner = 'TRIOmprage_1'
    # trg_seg_img_input_type = 'aseg'
    # trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([0, 2, 8, 9]))
    # trg_seg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_seg_img_input_type, np.arange(0,13))
    # all_seg_labels = []
    # for seg_file in trg_seg_filenames:
    #     print seg_file
    #     seg_img = nib.load(seg_file)
    #     seg_img_data = seg_img.get_data()
    #     seg_img_data_flatten = seg_img_data.flatten()
    #     seg_labels = np.unique(seg_img_data_flatten)
    #     all_seg_labels.append(seg_labels)
    #
    #
    # unique_seg_labels = np.unique(np.concatenate(all_seg_labels))
    # unique_seg_labels = unique_seg_labels[1:]
    # np.savetxt('aseg_labels.txt', unique_seg_labels, fmt='%d')
    #
    # trg_seg_img_input_type = 'aparc+aseg'
    # trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([0, 2, 8, 9]))
    # trg_seg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_seg_img_input_type, np.arange(0,13))
    # all_seg_labels = []
    # for seg_file in trg_seg_filenames:
    #     print seg_file
    #     seg_img = nib.load(seg_file)
    #     seg_img_data = seg_img.get_data()
    #     seg_img_data_flatten = seg_img_data.flatten()
    #     seg_labels = np.unique(seg_img_data_flatten)
    #     all_seg_labels.append(seg_labels)
    #
    #
    # unique_seg_labels = np.unique(np.concatenate(all_seg_labels))
    # unique_seg_labels = unique_seg_labels[1:]
    # np.savetxt('aparc+aseg_labels.txt', unique_seg_labels, fmt='%d')








