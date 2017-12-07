import os
import glob
import numpy as np
import tables
import nibabel as nib
import tempfile
from random import shuffle
from unet_model import  unet_model_3d
import sys
from os.path import join as opj

import subprocess
from keras import backend
from keras.models import load_model
from keras.optimizers import serialize
from keras.callbacks import Callback
from image_utils import patch_utils, intensity_standardize_utils
from nipype.interfaces.base import split_filename


class DeepImageSynthCallback(Callback):
    def __init__(self, output_prefix, save_per_epoch, save_weights, initial_epoch=1):
        self.output_prefix = output_prefix
        self.save_per_epoch = save_per_epoch
        self.save_weights = save_weights
        self.initial_epoch = initial_epoch
        super(DeepImageSynthCallback, self).__init__()

    def on_epoch_end(self, epoch, logs=None):
        sys.stdout.flush()
        print('')
        current_epoch = epoch + self.initial_epoch
        print('End of epoch %d' % current_epoch)
        if self.save_weights:
            weights_file = self.output_prefix + '_weights.h5'
            if self.save_per_epoch:
                root, ext = os.path.splitext(weights_file)
                weights_file = root + ('_epoch%d' % current_epoch) + ext
            self.model.save_weights(weights_file)
            print('Saving weights for epoch %d:' % current_epoch, weights_file)
        model_file = self.output_prefix + '_model.h5'
        if self.save_per_epoch:
            root, ext = os.path.splitext(model_file)
            model_file = root + ('_epoch%d' % current_epoch) + ext
        self.model.save(model_file)
        print('Saving model for epoch %d:' % current_epoch, model_file)
        print('')
        sys.stdout.flush()


def fetch_training_data_files(fs_dir, subjects_dir, img_input_type, training_subject_idxs):
    ''' assumes a freesurfer directory structure
    # Arguments
    :param fs_dir: directory with all the scanner freesurfer results stored
    :param src_scanner: scanner directory name from which we want to extract training images
    :param: src_img_input_type: freesurfer output that we want to extract for e.g. orig/001.mgz or aseg.mgz etc.
    :param trg_scanner
    :param trg_img_input_types tuple('orig/001', 'aparc+aseg') etc.
    '''
    src_subj_dir_list = sorted(glob.glob(os.path.join(fs_dir, subjects_dir, "[!fs]*", "mri")))
    training_data_files = list()
    input_subj_dir_list = list()
    for i in training_subject_idxs:
        input_subj_dir_list.append(src_subj_dir_list[i])
    for src_dir in input_subj_dir_list:
        training_data_files.append(os.path.join(src_dir, img_input_type + ".mgz"))
    return training_data_files


class DeepImageSynth(object):
    def __init__(self,  unet_num_filters, unet_depth, unet_downsampling_factor, feature_shape, storage_loc="memory", temp_folder=os.getcwd()):
        self.unet_downsampling_factor = unet_downsampling_factor
        self.unet_num_filters = unet_num_filters
        self.unet_depth = unet_depth
        self.feature_shape = feature_shape
        self.storage_loc = storage_loc
        self.temp_folder = temp_folder
        self.model = unet_model_3d(feature_shape, num_filters=unet_num_filters, unet_depth=unet_depth, downsize_filters_factor=unet_downsampling_factor,
                                              pool_size=(2, 2, 2), n_labels=1,initial_learning_rate=0.00001,
                                              deconvolution=False)

        self.model_trained = False
        self.model_compiled = False
        self.feature_generator = FeatureGenerator(self.feature_shape, temp_folder=temp_folder, storage_loc=storage_loc)
        self._weight_file = None

    @classmethod
    def from_file(cls, model_filename, storage_loc='memory', temp_folder=os.getcwd()):
        model = load_model(model_filename)
        input_shape = model.input_shape
        layer_list = model.get_config()["layers"]
        unet_num_filters = layer_list[1]['config']['filters']
        unet_downsampling_factor = 32/unet_num_filters
        num_pools = 0
        for layer in layer_list:
            if layer['class_name'] == 'MaxPooling3D':
                num_pools = num_pools + 1

        unet_depth = 2*num_pools - 1




        feature_shape = tuple(input_shape[1:-1])
        cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                       unet_downsampling_factor=unet_downsampling_factor,
                       feature_shape=feature_shape, storage_loc=storage_loc, temp_folder=temp_folder)
        cls_init.model = model
        cls_init.model_trained = True
        cls_init.model_compiled = True
        print('Loaded model file: %s' % model_filename)
        return cls_init

    def load_input_to_synth_images(self, image_filenames):
        if self.model is None:
            raise RuntimeError('Model does not exist')
        print('Extracting features (load_input_to_synth_images).')
        self.feature_generator.create_data_storage()
        self.feature_generator.create_feature_array(image_filenames, seg_filenames=None,
                                                    array_name='input_to_synth', indices=None)

    def load_training_images(self, source_filenames, target_filenames,
                             source_seg_filenames=None, target_seg_filenames=None, step_size=None):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_data(source_filenames, target_filenames,
                                                              source_seg_filenames, target_seg_filenames,
                                                              step_size=step_size)

    def load_validation_images(self, source_filenames, target_filenames,
                             source_seg_filenames=None, target_seg_filenames=None, step_size=None):
        if self.model is None:
            raise RuntimeError('Model does not exist')


        self.feature_generator.generate_src_trg_validation_data(source_filenames, target_filenames,
                                                              source_seg_filenames, target_seg_filenames,
                                                                step_size=step_size)

    def train_network(self, output_prefix, epochs=5, initial_epoch=1, batch_size=64, steps_per_epoch=10000,
                    optimizer='adam', loss='mean_squared_error', save_per_epoch=False, save_weights=True):
        print('Beginning Training. Using %s backend with "%s" data format.' % (backend.backend(),
                                                                               backend.image_data_format()))
        self.model.compile(optimizer, loss)
        if self._weight_file is not None:
            self.model.load_weights(self._weight_file)

        callback = DeepImageSynthCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)

        print('Training model...')
        self.model.fit_generator(generator=self.feature_generator.training_generator(batch_size=batch_size),
                                 epochs=epochs,
                                 validation_data=self.feature_generator.validation_generator(batch_size=batch_size),
                                 validation_steps=1000,
                                 steps_per_epoch=steps_per_epoch,
                                 callbacks=[callback], verbose=1, max_queue_size=100)
        self.model_trained = True

    def synthesize_image(self, in_img_file, out_img_filename, step_size):
        in_img = nib.load(in_img_file)
        (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
                                                                         step_size=step_size, is_label_img=False,
                                                                         indices=None)
        out_patches = self.model.predict(in_patches)
        patch_crop_size = [1, 1, 1]  # should be filter_size/2

        out_img_data, count_img_data = self.feature_generator.build_image_from_patches(out_patches, in_indices,
                                                                     padded_img_size, patch_crop_size)
        out_img_data = intensity_standardize_utils.wm_peak_normalize(out_img_data)
        out_img_data[out_img_data > 255] = 255

        out_img = nib.MGHImage(out_img_data, in_img.affine, in_img.header)
        nib.save(out_img, out_img_filename)


class FeatureGenerator(object):
    def __init__(self, feature_shape, temp_folder, storage_loc):
        self.feature_shape = feature_shape

        self.temp_folder = temp_folder
        self.storage_loc = storage_loc

        self.data_storage = None
        self.storage_filename = None

    def create_data_storage(self):
        if self.storage_loc == 'memory':
            self.storage = tables.open_file('tmp_data_storage.h5', 'w', driver='H5FD_CORE', driver_core_backing_store=False)
        elif self.storage_loc == 'disk':
            tmp_fp = tempfile.NamedTemporaryFile('w', suffix='.h5', dir=self.temp_folder, delete=False)
            self.storage_filename = tmp_fp.name
            tmp_fp.close()
            self.data_storage = tables.open_file(self.storage_filename, 'w')
        else:
            raise RuntimeError('Choose one of {memory, disk} for storage_loc')

    def load_data_storage(self, filename):
        if self.storage_loc== 'memory':
            self.data_storage = tables.open_file(filename, 'r', driver='core', driver_core_backing_store=False)
        elif self.storage_loc == 'disk':
            self.data_storage = tables.open_file(filename, 'r')
        else:
            raise RuntimeError('Choose one of {memory, disk} for storage_loc')

    def create_training_feature_array(self, image_filenames, seg_filenames, array_name, indices_list, step_size):

        nb_features_per_subject = 1000000
        nb_subjects = len(image_filenames)
        nb_src_modalities = len(image_filenames[0])
        print(image_filenames[0])
        tmpimg = nib.load(image_filenames[0])
        tmpseg = nib.load(seg_filenames[0])
        image_dtype = tmpimg.get_data_dtype()
        seg_dtype = tmpseg.get_data_dtype()

        feature_array = self.data_storage.create_earray(self.data_storage.root, array_name,
                                                 tables.Atom.from_dtype(image_dtype),
                                                 shape=(0,) + self.feature_shape + (1,),
                                                 expectedrows=np.prod(nb_features_per_subject)*nb_subjects)
        seg_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_seg',
                                                 tables.Atom.from_dtype(seg_dtype),
                                                 shape=(0,) + self.feature_shape + (1,),
                                                 expectedrows=np.prod(nb_features_per_subject)*nb_subjects)

        index_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_index',
                                                 tables.Int16Atom(), shape=(0, 3),
                                                 expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
        if indices_list == None:
            print("No indices_list found")
            indices_list = list()
            for input_file, seg_file in zip(image_filenames, seg_filenames):
                (features, indices) = self.extract_training_patches(input_file, seg_file, intensity_threshold=0,
                                                           step_size=step_size, indices=None)
                feature_array.append(features)
                index_array.append(indices)
                indices_list.append(indices)
                print(input_file + " features extract size ")
                print(features.shape)

        else:
            print("YES indices_list found")

            for input_file, seg_file, curr_indices in zip(image_filenames, seg_filenames, indices_list):
                print("curr indices shape is ")
                print(curr_indices.shape)
                (features, indices) = self.extract_training_patches(input_file, seg_file, intensity_threshold=0,
                                                           step_size=step_size, indices=curr_indices)

                print("indices shape is ")
                print(indices.shape)
                feature_array.append(features)
                index_array.append(curr_indices)
                print(input_file + " features extract size ")
                print(features.shape)



        return feature_array, index_array, indices_list

    def create_feature_array(self, image_filenames, seg_filenames, array_name, indices):

        nb_features_per_subject = 1000000
        nb_subjects = len(image_filenames)
        nb_src_modalities = len(image_filenames[0])
        print(image_filenames[0])
        tmpimg = nib.load(image_filenames[0])
        tmpseg = nib.load(seg_filenames[0])
        image_dtype = tmpimg.get_data_dtype()
        seg_dtype = tmpseg.get_data_dtype()

        feature_array = self.data_storage.create_earray(self.data_storage.root, array_name,
                                                 tables.Atom.from_dtype(image_dtype),
                                                 shape=(0,) + self.feature_shape + (1,),
                                                 expectedrows=np.prod(nb_features_per_subject)*nb_subjects)
        seg_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_seg',
                                                 tables.Atom.from_dtype(seg_dtype),
                                                 shape=(0,) + self.feature_shape + (1,),
                                                 expectedrows=np.prod(nb_features_per_subject)*nb_subjects)

        index_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_index',
                                                 tables.Int16Atom(), shape=(0, 3),
                                                 expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        for input_file, seg_file in zip(image_filenames, seg_filenames):
            if indices == None:
                (features, indices) = self.extract_patches(input_file, intensity_threshold=0,
                                                           step_size=[1,1,1], indices=None)

            else:
                (features, indices) = self.extract_patches(input_file, intensity_threshold=0,
                                                           step_size=[1,1,1], indices=indices)
            print(features.shape)
            print(indices.shape)

            print(feature_array.shape)
            print(index_array.shape)

            feature_array.append(features)
            index_array.append(indices)

        return feature_array, index_array, indices

    def extract_training_patches(self, in_img_file, seg_img_file, intensity_threshold, step_size, indices):

        if indices is not None:
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                  is_label_img=False, indices=indices)


            return patches, indices
        else:
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                      is_label_img=False, indices=indices)
            (seg_patches, seg_indices, _) = self.extract_patches(seg_img_file, intensity_threshold, step_size,
                                                              is_label_img=True, indices=indices)
            training_patches = patches
            training_indices = indices
            # print("collecting patches from " + in_img_file)
            # #training_patches = np.zeros((800000, self.feature_shape[0], self.feature_shape[0], self.feature_shape[0], 1))
            # training_patches = []
            #
            # print("collecting patches from " + seg_img_file)
            # #training_indices = np.zeros((800000, 3))
            # training_indices = []
            # pointer = 0
            # center_voxel_seg = seg_patches[:,self.feature_shape[0]/2,self.feature_shape[1]/2,self.feature_shape[2]/2,0]
            # unique_labels = np.unique(center_voxel_seg)
            # min_samples_per_label =  5000
            # for l in unique_labels:
            #     print(l)
            #     if l == 0:
            #         continue
            #     elif l == 7 | l == 8 | l == 15 | l == 16 | l == 46 | l == 47:
            #         # cerebellum. Choose all voxels in training
            #         l_idxs = np.where(center_voxel_seg == l)[0]
            #         # choose 10% of total number of label voxels or 5000
            #         percent_samples = len(l_idxs) / 1
            #         samples_per_label = np.maximum(percent_samples, 10000)
            #         rand_idx_idxs = np.random.randint(0, len(l_idxs), (samples_per_label))
            #         chosen_patches = patches[l_idxs[rand_idx_idxs]]
            #         chosen_indices = indices[l_idxs[rand_idx_idxs]]
            #         #training_patches[pointer:pointer + len(chosen_patches)] = chosen_patches
            #         training_patches.append(chosen_patches)
            #         #training_indices[pointer:pointer + len(chosen_patches)] = chosen_indices
            #         training_indices.append(chosen_indices)
            #         pointer = pointer + len(chosen_patches)
            #
            #     else:
            #         l_idxs = np.where(center_voxel_seg == l)[0]
            #         # choose 10% of total number of label voxels or 5000
            #         onepercent_samples = len(l_idxs) / 10
            #         samples_per_label = np.maximum(onepercent_samples, min_samples_per_label)
            #         rand_idx_idxs = np.random.randint(0, len(l_idxs), (samples_per_label))
            #         chosen_patches = patches[l_idxs[rand_idx_idxs]]
            #         chosen_indices = indices[l_idxs[rand_idx_idxs]]
            #         #training_patches[pointer:pointer + len(chosen_patches)] = chosen_patches
            #         training_patches.append(chosen_patches)
            #         #training_indices[pointer:pointer + len(chosen_patches)] = chosen_indices
            #         training_indices.append(chosen_indices)
            #
            #         pointer = pointer + len(chosen_patches)
            #
            #
            # #training_patches = training_patches[0:pointer]
            # #training_indices = training_indices[0:pointer]
            # training_patches = np.array(training_patches)
            # training_indices = np.array(training_indices)

            return training_patches, np.int32(training_indices)


    def extract_patches(self, in_img_file, intensity_threshold, step_size, is_label_img=False, indices=None):
        # pad the images by patch shape

        in_img = nib.load(in_img_file)
        in_img_data = in_img.get_data()

        # white matter peak set to 200 and divide by 255
        if is_label_img == False:
            in_img_data = intensity_standardize_utils.wm_peak_normalize(in_img_data)
            in_img_data[in_img_data > 255] = 255
            in_img_data = in_img_data / 255




        in_img_data_pad = np.pad(in_img_data, ((self.feature_shape[0], self.feature_shape[0]),
                                               (self.feature_shape[1], self.feature_shape[1]),
                                               (self.feature_shape[2], self.feature_shape[2])),
                                 'constant', constant_values=0)

        padded_img_size = in_img_data_pad.shape

        if  indices is not None :

            idx_x = indices[:, 0]
            idx_y = indices[:, 1]
            idx_z = indices[:, 2]
        else:
            (idx_x_fg, idx_y_fg, idx_z_fg) = np.where(in_img_data_pad > intensity_threshold)
            min_idx_x_fg = np.min(idx_x_fg) - step_size[0]
            max_idx_x_fg = np.max(idx_x_fg) + step_size[0]
            min_idx_y_fg = np.min(idx_y_fg) - step_size[1]
            max_idx_y_fg = np.max(idx_y_fg) + step_size[1]
            min_idx_z_fg = np.min(idx_z_fg) - step_size[2]
            max_idx_z_fg = np.max(idx_z_fg) + step_size[2]

            sampled_x = np.arange(min_idx_x_fg, max_idx_x_fg, step_size[0])
            sampled_y = np.arange(min_idx_y_fg, max_idx_y_fg, step_size[1])
            sampled_z = np.arange(min_idx_z_fg, max_idx_z_fg, step_size[2])

            idx_x, idx_y, idx_z = np.meshgrid(sampled_x, sampled_y, sampled_z, sparse=False, indexing='ij')
            idx_x = idx_x.flatten()
            idx_y = idx_y.flatten()
            idx_z = idx_z.flatten()

        patches = []

        for patch_iter in range(len(idx_x)):
            curr_patch = in_img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                         idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1],
                         idx_z[patch_iter]:idx_z[patch_iter] + self.feature_shape[2]]

            patches.append(curr_patch)

        patches = np.asarray(patches)

        # add channel as a dimension for keras
        newshape = list(patches.shape)
        newshape.append(1)
        # print newshape
        # print newshape.__class__
        patches = np.reshape(patches, newshape)
        indices = np.concatenate((idx_x.reshape(-1,1), idx_y.reshape(-1, 1), idx_z.reshape(-1, 1)), axis=1)
        return patches, np.int32(indices), padded_img_size

    def build_image_from_patches(self, in_patches, indices, padded_img_size, patch_crop_size):
        ''' patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]'''
        import numpy as np
        out_img_data = np.zeros(padded_img_size)
        count_img_data = np.zeros(padded_img_size)
        patch_mask = np.zeros(self.feature_shape)
        patch_mask[0 + patch_crop_size[0]: self.feature_shape[0] - patch_crop_size[0],
        0 + patch_crop_size[1]: self.feature_shape[1] - patch_crop_size[1],
        0 + patch_crop_size[2]: self.feature_shape[2] - patch_crop_size[2]] = 1

        idx_x = indices[:,0]
        idx_y = indices[:,1]
        idx_z = indices[:,2]

        for patch_iter in range(len(idx_x)):
            out_img_data[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + self.feature_shape[2]] += \
                np.multiply(np.reshape(in_patches[patch_iter, :],in_patches.shape[1:4]), patch_mask)

            count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + self.feature_shape[2]] += patch_mask

        out_img_data = np.divide(out_img_data, count_img_data)
        out_img_data[np.isnan(out_img_data)] = 0
        # remove the padding
        unpadded_img_size = padded_img_size - np.multiply(self.feature_shape, 2)

        out_img_data = out_img_data[self.feature_shape[0]:self.feature_shape[0] + unpadded_img_size[0],
                       self.feature_shape[1]:self.feature_shape[1] + unpadded_img_size[1],
                       self.feature_shape[2]:self.feature_shape[2] + unpadded_img_size[2]]
        count_img_data = count_img_data[self.feature_shape[0]:self.feature_shape[0] + unpadded_img_size[0],
                         self.feature_shape[1]:self.feature_shape[1] + unpadded_img_size[1],
                         self.feature_shape[2]:self.feature_shape[2] + unpadded_img_size[2]]

        return out_img_data, count_img_data


    def generate_src_trg_training_data(self, source_filenames, target_filenames,
                                       source_seg_filenames=None, target_seg_filenames=None, step_size=None):
        print('Creating source image patches.')
        (_,_, indices_list) = self.create_training_feature_array(source_filenames, source_seg_filenames, 'src',
                                                                 indices_list=None, step_size=step_size)

        print('Creating target image patches.')
        self.create_training_feature_array(target_filenames, target_seg_filenames, 'trg', indices_list,
                                           step_size=step_size)

    def generate_src_trg_validation_data(self, source_filenames, target_filenames,
                                       source_seg_filenames=None, target_seg_filenames=None, step_size=None):
        print('Creating source image patches.')
        (_,_, indices_list) = self.create_training_feature_array(source_filenames, source_seg_filenames,
                                                            'src_validation', indices_list=None, step_size=step_size)

        print('Creating target image patches.')
        self.create_training_feature_array(target_filenames, target_seg_filenames,
                                           'trg_validation', indices_list=indices_list, step_size=step_size)

    def training_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src[index])
                y_list.append(self.data_storage.root.trg[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list

    def validation_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src_validation.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src_validation[index])
                y_list.append(self.data_storage.root.src_validation[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list




if __name__ == "__main__":


    fs_dir = '/autofs/space/mreuter/users/amod/pfizer_dataset_analysis/data/fs_syn_reg_dir_v1/freesurfer6p0_skullstripped_v1'
    src_img_input_type = 'orig/001'
    src_scanner = 'TRIOmecho'
    src_filenames = fetch_training_data_files(fs_dir, src_scanner, src_img_input_type, np.array([0,2,8,9]))
    src_seg_img_input_type = 'aparc+aseg'
    src_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([0,2,8,9]))

    trg_img_input_type = 'orig/001'
    trg_scanner = 'TRIOmprage_1'
    trg_seg_img_input_type = 'aparc+aseg'
    trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([0,2,8,9]))
    trg_seg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_seg_img_input_type, np.array([0,2,8,9]))

    src_img_input_type = 'orig/001'
    src_scanner = 'TRIOmecho'
    val_src_filenames = fetch_training_data_files(fs_dir, src_scanner, src_img_input_type, np.array([3]))
    src_seg_img_input_type = 'aparc+aseg'
    val_src_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([3]))

    trg_img_input_type = 'orig/001'
    trg_scanner = 'TRIOmprage_1'
    trg_seg_img_input_type = 'aparc+aseg'
    val_trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([3]))
    val_trg_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([3]))

    fg = FeatureGenerator(feature_shape = (64, 64, 64),
                          temp_folder="/local_mount/space/nike/1/users/amod/tmp", storage_loc="disk")
    fg.create_data_storage()
    fg.generate_src_trg_training_data(src_filenames, trg_filenames, src_seg_filenames, trg_seg_filenames)
    fg.generate_src_trg_validation_data(val_src_filenames, val_trg_filenames, val_src_seg_filenames, val_trg_seg_filenames)
