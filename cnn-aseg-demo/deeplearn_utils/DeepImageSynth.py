import os
import glob
import numpy as np
import tables
import nibabel as nib
import tempfile
from random import shuffle
from unet_model import  unet_model_3d, atrous_net, grad_loss, noise_net, class_net, \
    pure_grad_loss, dice_coef_loss, dice_coef_loss2
import sys
from os.path import join as opj

import subprocess
from keras import backend
from keras.models import load_model
from keras.optimizers import serialize
from keras.callbacks import Callback
from image_utils.image_utils import patch_utils, intensity_standardize_utils
from keras.models import  Model
from keras.callbacks import ReduceLROnPlateau, TensorBoard, ModelCheckpoint
from sklearn import preprocessing

def detach_model(m):
    for l in m.layers:
        if l.name == 'model_1':
            return l
    return m


class MultiGPUCheckpointCallback(Callback):

    def __init__(self, output_prefix, save_per_epoch, save_weights, initial_epoch=1):
        super(MultiGPUCheckpointCallback, self).__init__()
        self.output_prefix = output_prefix
        self.save_per_epoch = save_per_epoch
        self.save_weights = save_weights
        self.initial_epoch = initial_epoch
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

            detach_model(self.model).save_weights(weights_file)
            print('Saving weights for epoch %d:' % current_epoch, weights_file)
        model_file = self.output_prefix + '_model.h5'
        if self.save_per_epoch:
            root, ext = os.path.splitext(model_file)
            model_file = root + ('_epoch%d' % current_epoch) + ext
        detach_model(self.model).save(model_file)
        print('Saving model for epoch %d:' % current_epoch, model_file)
        print('')
        sys.stdout.flush()

    # def on_epoch_end(self, epoch, logs=None):
    #     logs = logs or {}
    #     self.epochs_since_last_save += 1
    #     if self.epochs_since_last_save >= self.period:
    #         self.epochs_since_last_save = 0
    #         filepath = self.filepath.format(epoch=epoch + 1, **logs)
    #         if self.verbose > 0:
    #             print('Epoch %05d: saving model to %s' % (epoch + 1, filepath))
    #             self.base_model.save_weights(filepath, overwrite=True)
    #             self.base_model.save(filepath, overwrite=True)

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
    def __init__(self,  unet_num_filters, unet_depth, unet_downsampling_factor, feature_shape, dim=3,
                 storage_loc="memory", temp_folder=os.getcwd(), n_labels=0, labels=[], net='unet', loss='mean_absolute_error',
                 initial_learning_rate=0.00001, use_patches=True, wmp_standardize=True, fcn=True, num_gpus=1, preprocessing=False):
        self.net = net
        self.feature_shape = feature_shape
        self.storage_loc = storage_loc
        self.temp_folder = temp_folder
        self.wmp_standardize = wmp_standardize
        self.fcn = fcn # fully convolutional or not boolean
        self.num_gpus = num_gpus
        self.preprocessing = preprocessing


        if net == 'unet':
            self.unet_downsampling_factor = unet_downsampling_factor
            self.unet_num_filters = unet_num_filters
            self.unet_depth = unet_depth

            self.model, self.parallel_model = unet_model_3d(feature_shape, num_filters=unet_num_filters, unet_depth=unet_depth,
                                       downsize_filters_factor=unet_downsampling_factor, pool_size=(2, 2, 2),
                                       n_labels=n_labels,loss=loss, initial_learning_rate=initial_learning_rate,
                                       deconvolution=False, use_patches=use_patches, num_gpus=num_gpus)
        elif net == 'atrousnet':
            self.model = atrous_net(feature_shape, unet_num_filters, initial_learning_rate=initial_learning_rate,
                                    loss='mean_absolute_error')
        elif net == 'noisenet':
            self.model = noise_net(feature_shape, unet_num_filters, initial_learning_rate=initial_learning_rate,
                                    loss='mean_absolute_error')
        elif (net == 'class_net') & (fcn==False):
            # i.e. not a FCN like the above. A single label for the input feature, not a patch of labels as above
            self.model = class_net(feature_shape=feature_shape, dim=dim, unet_num_filters=unet_num_filters,
                                   n_labels=n_labels, initial_learning_rate=initial_learning_rate, loss=loss)




        self.model_trained = False
        self.model_compiled = False
        self.labels = labels
        self.n_labels = len(labels)
        self.feature_generator = FeatureGenerator(self.feature_shape, temp_folder=temp_folder,
                                                  storage_loc=storage_loc, labels=labels, n_labels=n_labels,
                                                  wmp_standardize=wmp_standardize, use_patches=use_patches, dim=dim,
                                                  preprocessing=self.preprocessing)
        self._weight_file = None

    @classmethod
    def from_file(cls, model_filename, loss, storage_loc='memory', temp_folder=os.getcwd(),
                  net='unet', n_labels=0, labels=None):
        print loss
        if loss == 'dice_coef_loss':
            model = load_model(model_filename,custom_objects={'dice_coef_loss': dice_coef_loss})
        elif loss == 'dice_coef_loss2':
            model = load_model(model_filename, custom_objects={'dice_coef_loss2': dice_coef_loss2})
        elif loss == 'grad_loss':
            model = load_model(model_filename, custom_objects={'grad_loss': grad_loss})
        else:
            model = load_model(model_filename)

        input_shape = model.input_shape
        layer_list = model.get_config()["layers"]
        unet_num_filters = layer_list[1]['config']['filters']
        unet_downsampling_factor = unet_num_filters / unet_num_filters
        unet_loss = model.loss
        dim = len(input_shape) - 2

        num_pools = 0
        for layer in layer_list:
            if layer['class_name'] == 'MaxPooling3D':
                num_pools = num_pools + 1

        unet_depth = 2 * num_pools - 1

        feature_shape = tuple(input_shape[1:-1])

        if net=='unet':

            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                       unet_downsampling_factor=unet_downsampling_factor,
                       feature_shape=feature_shape, loss=loss, storage_loc=storage_loc,
                        n_labels=n_labels,labels=labels,
                        temp_folder=temp_folder, net='unet')
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True
            print('Loaded model file: %s' % model_filename)
        elif net=='atrousnet':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                       unet_downsampling_factor=unet_downsampling_factor,
                       feature_shape=feature_shape, storage_loc=storage_loc, temp_folder=temp_folder, net='atrousnet')
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True
            print('Loaded model file: %s' % model_filename)
        elif net == 'noisenet':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                           unet_downsampling_factor=unet_downsampling_factor,
                           feature_shape=feature_shape, storage_loc=storage_loc, temp_folder=temp_folder,
                           net='noisenet')
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True
            print('Loaded model file: %s' % model_filename)


        return cls_init

    @classmethod
    def from_file_old_model(cls, model_filename, loss, storage_loc='memory', temp_folder=os.getcwd(),
                  net='class_net', n_labels=0):
        model = load_model(model_filename)

        input_shape = model.input_shape
        unet_num_filters = model.get_config()[0]['config']['filters']
        dim = len(input_shape) - 2



        unet_depth = len(model.get_config()[0])/4

        feature_shape = tuple(input_shape[1:-1])
        if net == 'class_net':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth, unet_downsampling_factor=1,
                           feature_shape=feature_shape, dim=dim,storage_loc=storage_loc,
                           temp_folder=temp_folder, n_labels=n_labels,
                           net='class_net')
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True

        return cls_init




    # def load_input_to_synth_images(self, image_filenames, is_label_img):
    #     if self.model is None:
    #         raise RuntimeError('Model does not exist')
    #     print('Extracting features (load_input_to_synth_images).')
    #     self.feature_generator.create_data_storage()
    #     self.feature_generator.create_feature_array(image_filenames, array_name='input_to_synth', indices=None, is_label_img=is_label_img)

    def load_training_images(self, source_filenames, target_filenames,is_src_label_img, is_trg_label_img,
                             source_seg_filenames=None, target_seg_filenames=None, step_size=None,
                             preprocessing=False
                             ):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_data(source_filenames, target_filenames,
                                                              is_src_label_img, is_trg_label_img,
                                                              step_size=step_size,
                                                              )

    def load_training_slices_and_labels(self, source_filenames, label_list, is_src_label_img):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_data(source_filenames, target_filenames=None, is_trg_label_img=False,
                                                              is_src_label_img=is_src_label_img, target_label_list=label_list)




    def load_validation_images(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                             source_seg_filenames=None, target_seg_filenames=None, step_size=None,preprocessing=False):
        if self.model is None:
            raise RuntimeError('Model does not exist')


        self.feature_generator.generate_src_trg_validation_data(source_filenames, target_filenames,
                                                                is_src_label_img, is_trg_label_img,
                                                                step_size=step_size)


    def load_validation_slices_and_labels(self, source_filenames, label_list, is_src_label_img):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.generate_src_trg_validation_data(source_filenames, target_filenames=None, is_trg_label_img=False,
                                                              is_src_label_img=is_src_label_img, target_label_list=label_list)



    def train_network(self, output_prefix, epochs=5, initial_epoch=1, batch_size=64, steps_per_epoch=10000, validation_steps=1000,
                    optimizer='adam', save_per_epoch=False, save_weights=True):
        print('Beginning Training. Using %s backend with "%s" data format.' % (backend.backend(),
                                                                               backend.image_data_format()))
        if self._weight_file is not None:
            self.model.load_weights(self._weight_file)

        reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5 ,patience=5, min_lr=0.000001)
        # modelcp = ModelCheckpoint(output_prefix, monitor='val_loss', verbose=0, save_best_only=False,
        #                                 save_weights_only=False, mode='auto', period=1)

        print('Training model...')
        if self.n_labels == 0:
            if self.num_gpus == 1:
                callback = DeepImageSynthCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)

                self.model.fit_generator(generator=self.feature_generator.training_generator(batch_size=batch_size*self.num_gpus),
                                         epochs=epochs,
                                         validation_data=self.feature_generator.validation_generator(batch_size=batch_size*self.num_gpus),
                                         validation_steps=validation_steps,
                                         steps_per_epoch=steps_per_epoch,
                                         callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)
                self.model_trained = True
            else:
                callback = MultiGPUCheckpointCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)

                self.parallel_model.fit_generator(
                    generator=self.feature_generator.training_generator(batch_size=batch_size * self.num_gpus),
                    epochs=epochs,
                    validation_data=self.feature_generator.validation_generator(batch_size=batch_size * self.num_gpus),
                    validation_steps=validation_steps,
                    steps_per_epoch=steps_per_epoch,
                    callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=1000)
                self.parallel_model_trained = True

        elif (self.n_labels > 1) & (self.fcn == True) :
            if self.num_gpus == 1:
                callback = DeepImageSynthCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)
                self.model.fit_generator(generator=self.feature_generator.seg_training_generator(batch_size=batch_size*self.num_gpus),
                                     epochs=epochs,
                                     validation_data=self.feature_generator.seg_validation_generator(batch_size=batch_size*self.num_gpus),
                                     validation_steps=validation_steps,
                                     steps_per_epoch=steps_per_epoch,
                                     callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=10)
            else:
                callback = MultiGPUCheckpointCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)
                self.parallel_model.fit_generator(
                    generator=self.feature_generator.training_generator(batch_size=batch_size * self.num_gpus),
                    epochs=epochs,
                    validation_data=self.feature_generator.validation_generator(batch_size=batch_size * self.num_gpus),
                    validation_steps=validation_steps,
                    steps_per_epoch=steps_per_epoch,
                    callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)
                self.parallel_model_trained = True



        elif (self.n_labels > 1) & (self.fcn == False):
            self.model.fit_generator(generator=self.feature_generator.training_label_generator(batch_size=batch_size*self.num_gpus),
                                     epochs=epochs,
                                     validation_data=self.feature_generator.validation_label_generator(
                                         batch_size=batch_size * self.num_gpus),
                                     validation_steps=validation_steps,
                                     steps_per_epoch=steps_per_epoch,
                                     callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)



    def synthesize_image(self, in_img_file, out_img_filename, step_size):
        in_img = nib.load(in_img_file)
        print("Shape is " + str(in_img.get_data().shape))
        (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
                                                                         step_size=step_size, is_label_img=False,
                                                                         indices=None)

        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(in_patches)
        else:
            out_patches = self.model.predict(in_patches)

        patch_crop_size = [1, 1, 1]  # should be filter_size/2

        print("padded image size " + str(padded_img_size))
        out_img_data, count_img_data = self.feature_generator.build_image_from_patches(out_patches, in_indices,
                                                                     padded_img_size, patch_crop_size, step_size)
        out_img_data = intensity_standardize_utils.wm_peak_normalize(out_img_data)
        print("Out data shape is " + str(out_img_data.shape))
        # out_img_data[out_img_data > 255] = 255

        out_img = nib.MGHImage(out_img_data, in_img.affine, in_img.header)
        nib.save(out_img, out_img_filename)

    def predict_segmentation(self, in_img_file, out_soft_filename, out_hard_filename, step_size):
        in_img = nib.load(in_img_file)
        (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
                                                                         step_size=step_size, is_label_img=False,
                                                                         indices=None)

        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(in_patches)
        else:
            out_patches = self.model.predict(in_patches)



        num_labels = out_patches.shape[-1]
        padded_img_size_multiple_labels = padded_img_size + (num_labels,)


        patch_crop_size = [1, 1, 1]  # should be filter_size/2

        out_img_data, label_img_data, count_img_data = self.feature_generator.build_seg_from_patches(out_patches, in_indices,
                                                                     padded_img_size_multiple_labels, patch_crop_size, step_size)



        out_img = nib.MGHImage(out_img_data, in_img.affine, in_img.header)
        nib.save(out_img, out_soft_filename)

        label_img = nib.MGHImage(label_img_data, in_img.affine, in_img.header)
        nib.save(label_img, out_hard_filename)

    def predict_labels(self, in_img_file):
        in_img = nib.load(in_img_file)
        (in_slices, in_indices) = self.feature_generator.extract_slices(in_img_file, intensity_threshold=0,
                                                                        is_label_img=False,indices=None)

        out_prob = self.model.predict(in_slices)
        out_labels = np.zeros(out_prob.shape)
        out_labels[out_prob > 0.5] = 1

        return out_prob, out_labels

    def apply_encoder(self, in_img_file, layer_name, step_size):
        encoder_model = Model(inputs=self.model.input, outputs=self.model.get_layer(layer_name).output)
        in_img = nib.load(in_img_file)
        print("Shape is " + str(in_img.get_data().shape))
        (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
                                                                         step_size=step_size, is_label_img=False,
                                                                         indices=None)

        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(in_patches)
        else:
            out_patches = encoder_model.predict(in_patches)

        return out_patches, in_patches,  in_indices, padded_img_size

    def apply_decoder(self, in_patches, layer_name, step_size):
        decoder_model = Model(inputs=self.model.get_layer(layer_name), outputs=self.model.output)
        # in_img = nib.load(in_img_file)
        # print("Shape is " + str(in_img.get_data().shape))
        # (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
        #                                                                  step_size=step_size, is_label_img=False,
        #                                                                  indices=None)

        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(in_patches)
        else:
            out_patches = decoder_model.predict(in_patches)

        return out_patches



class FeatureGenerator(object):
    def __init__(self, feature_shape, temp_folder, storage_loc, labels=None, n_labels=0,
                 wmp_standardize=True, use_patches=False, dim=3, preprocessing=False):

        self.feature_shape = feature_shape

        self.temp_folder = temp_folder
        self.storage_loc = storage_loc

        self.n_labels = n_labels
        self.labels = labels
        self.wmp_standardize = wmp_standardize
        self.use_patches = use_patches
        self.dim = dim
        self.preprocessing = preprocessing


        self.data_storage = None
        self.storage_filename = None


    def create_data_storage(self):
        if self.storage_loc == 'memory':
            self.data_storage = tables.open_file('tmp_data_storage.h5', 'w', driver='H5FD_CORE', driver_core_backing_store=False)
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

    def create_training_feature_array(self, image_filenames,  array_name, indices_list, step_size, is_label_img):
        nb_features_per_subject = 1000000
        nb_subjects = len(image_filenames)
        nb_src_modalities = len(image_filenames[0])
        print(image_filenames[0])
        tmpimg = nib.load(image_filenames[0])
        # tmpseg = nib.load(seg_filenames[0])
        if is_label_img == True:
            image_dtype = tmpimg.get_data_dtype()
        else:
            image_dtype = np.dtype(np.float32)

        feature_array = self.data_storage.create_earray(self.data_storage.root, array_name,
                                                        tables.Atom.from_dtype(image_dtype),
                                                        shape=(0,) + self.feature_shape + (1,),
                                                        expectedrows=np.prod(nb_features_per_subject) * nb_subjects)



        if self.use_patches == True:
            index_array = self.data_storage.create_earray(self.data_storage.root, array_name + '_index',
                                                          tables.Int16Atom(), shape=(0, 3),
                                                          expectedrows=np.prod(nb_features_per_subject) * nb_subjects)


            if indices_list == None:
                print("No indices_list found")
                indices_list = list()
                for input_file in image_filenames:
                    (features, indices) = self.extract_training_patches(input_file, intensity_threshold=0,
                                                               step_size=step_size, indices=None, is_label_img=is_label_img)
                    feature_array.append(features)
                    index_array.append(indices)
                    indices_list.append(indices)
                    print(input_file + " features extract size ")
                    print(features.shape)

            else:
                print("YES indices_list found")

                for input_file, curr_indices in zip(image_filenames, indices_list):
                    print("curr indices shape is ")
                    print(curr_indices.shape)
                    (features, indices) = self.extract_training_patches(input_file, intensity_threshold=0,
                                                               step_size=step_size, indices=curr_indices,
                                                                        is_label_img=is_label_img)

                    print("indices shape is ")
                    print(indices.shape)
                    feature_array.append(features)
                    index_array.append(curr_indices)
                    print(input_file + " features extract size ")
                    print(features.shape)
        else:

            index_array = self.data_storage.create_earray(self.data_storage.root, array_name + '_index',
                                                          tables.Int16Atom(), shape=(0, 1),
                                                          expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
            if indices_list == None:
                indices_list = list()

                if self.dim == 2:
                    for input_file in image_filenames:
                        (features, indices) = self.extract_training_slices(input_file, intensity_threshold=0,
                                                                           indices=None,is_label_img=is_label_img)
                        feature_array.append(features)
                        index_array.append(indices)
                        indices_list.append(indices)
                        print(input_file + " features extract size ")
                        print(features.shape)
                else:
                    print("YES indices_list found")

                    for input_file, curr_indices in zip(image_filenames, indices_list):
                        print("curr indices shape is ")
                        print(curr_indices.shape)
                        (features, indices) = self.extract_training_slices(input_file, intensity_threshold=0,
                                                                           indices=curr_indices, is_label_img=is_label_img)

                        print("indices shape is ")
                        print(indices.shape)
                        feature_array.append(features)
                        index_array.append(curr_indices)
                        print(input_file + " features extract size ")
                        print(features.shape)


        return feature_array, index_array, indices_list

    def create_training_label_array(self, target_label_list,  array_name, indices_list):
        nb_features_per_subject = 200
        nb_subjects = len(target_label_list)


        label_array = self.data_storage.create_earray(self.data_storage.root, array_name,
                                                      tables.Int16Atom(), shape=(0, 1),
                                                        expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
        index_array = self.data_storage.create_earray(self.data_storage.root, array_name + '_index',
                                                      tables.Int16Atom(), shape=(0, 1),
                                                      expectedrows=np.prod(nb_features_per_subject) * nb_subjects)


        for curr_label, curr_indices in zip(target_label_list, indices_list):
            print("curr indices shape is ")
            print(curr_indices.shape)

            curr_label_array = np.ones(curr_indices.shape)*curr_label
            label_array.append(curr_label_array)
            index_array.append(curr_indices)


    # def create_feature_array(self, image_filenames, array_name, indices, is_label_img):
    #
    #     nb_features_per_subject = 1000000
    #     nb_subjects = len(image_filenames)
    #     nb_src_modalities = len(image_filenames[0])
    #     print(image_filenames[0])
    #     tmpimg = nib.load(image_filenames[0])
    #     # tmpseg = nib.load(seg_filenames[0])
    #     image_dtype = np.dtype(np.float32) #tmpimg.get_data_dtype()
    #     # seg_dtype = tmpseg.get_data_dtype()
    #
    #     feature_array = self.data_storage.create_earray(self.data_storage.root, array_name,
    #                                              tables.Atom.from_dtype(image_dtype),
    #                                              shape=(0,) + self.feature_shape + (1,),
    #                                              expectedrows=np.prod(nb_features_per_subject)*nb_subjects)
    #     # seg_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_seg',
    #     #                                          tables.Atom.from_dtype(seg_dtype),
    #     #                                          shape=(0,) + self.feature_shape + (1,),
    #     #                                          expectedrows=np.prod(nb_features_per_subject)*nb_subjects)
    #
    #     if self.use_patches == True:
    #         index_array = self.data_storage.create_earray(self.data_storage.root, array_name+'_index',
    #                                                  tables.Int16Atom(), shape=(0, 3),
    #                                                  expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
    #
    #         for input_file in image_filenames:
    #             if indices == None:
    #                 (features, indices,_) = self.extract_patches(input_file, intensity_threshold=0,is_label_img=is_label_img,
    #                                                            step_size=[1,1,1], indices=None)
    #
    #             else:
    #                 (features, indices, _) = self.extract_patches(input_file, intensity_threshold=0,
    #                                                            step_size=[1,1,1], indices=indices)
    #             print(features.shape)
    #             print(indices.shape)
    #
    #             print(feature_array.shape)
    #             print(index_array.shape)
    #
    #             feature_array.append(features)
    #             index_array.append(indices)
    #     else:
    #
    #         index_array = self.data_storage.create_earray(self.data_storage.root, array_name + '_index',
    #                                                       tables.Int16Atom(), shape=(0, 1),
    #                                                       expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
    #         if self.dim == 2:
    #             for input_file in image_filenames:
    #                 (features, indices) = self.extract_training_slices(input_file, intensity_threshold=0,
    #                                                                    indices=None, is_label_img=is_label_img)
    #                 feature_array.append(features)
    #                 index_array.append(indices)
    #                 indices_list.append(indices)
    #
    #     return feature_array, index_array, indices

    def extract_training_slices(self,in_img_file, intensity_threshold, indices, is_label_img ):

        if indices is not None:
            (slices, indices) = self.extract_slices(in_img_file, intensity_threshold,
                                                         is_label_img=is_label_img, indices=indices)

            return slices, indices
        else:
            (slices, indices) = self.extract_slices(in_img_file, intensity_threshold,
                                                    is_label_img=is_label_img, indices=indices)

            training_slices = slices
            training_indices = indices


            return training_slices, np.int32(training_indices)

    def extract_slices(self, in_img_file, intensity_threshold, indices=None, is_label_img=False ):
        in_img = nib.load(in_img_file)


        # white matter peak set to 200 and divide by 255
        if is_label_img == False:
            in_img_data = in_img.get_data().astype(float)

            if self.wmp_standardize == True:
                in_img_data = intensity_standardize_utils.wm_peak_normalize(in_img_data)
                in_img_data[in_img_data > 255] = 255
                in_img_data = in_img_data / 255
            else:
                in_img_data = intensity_standardize_utils.robust_normalize(in_img_data)
                in_img_data = in_img_data / 255


            print("Image max is :" + str(in_img_data.max()))
        else:
            in_img_data = in_img.get_data()

        if indices is not None:
            slices = []
            for z_index in indices:
                slices.append(in_img_data[:,:,z_index])


        else:
            slices = []
            indices = range(in_img_data.shape[2])
            print(in_img_data.shape[2])
            for z_index in indices :
                print(z_index)
                slices.append(in_img_data[:, :, z_index])

            slices = np.asarray(slices)
            indices = np.asarray(indices)

            # add channel as a dimension for keras
            newshape = list(slices.shape)
            newshape.append(1)
            # print newshape
            # print newshape.__class__
            slices = np.reshape(slices, newshape)
            indices = indices.reshape(-1, 1)

        return slices, indices







    def extract_training_patches(self, in_img_file, intensity_threshold, step_size, indices, is_label_img):

        if indices is not None:
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                  is_label_img=is_label_img, indices=indices)


            return patches, indices
        else:
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                      is_label_img=is_label_img, indices=indices)
            # (seg_patches, seg_indices, _) = self.extract_patches(seg_img_file, intensity_threshold, step_size,
            #                                                   is_label_img=True, indices=indices)
            training_patches = patches
            training_indices = indices


            return training_patches, np.int32(training_indices)


    def extract_patches(self, in_img_file, intensity_threshold, step_size, is_label_img=False, indices=None):
        # pad the images by patch shape

        in_img = nib.load(in_img_file)


        # white matter peak set to 200 and divide by 255
        if is_label_img == False:
            in_img_data = in_img.get_data().astype(float)


            if self.wmp_standardize == True:
                in_img_data = intensity_standardize_utils.wm_peak_normalize(in_img_data)
                in_img_data[in_img_data > 255] = 255
                in_img_data = in_img_data / 255


            else:
                in_img_data = intensity_standardize_utils.robust_normalize(in_img_data)
                in_img_data = in_img_data / 255

        else:
            in_img_data = in_img.get_data()
            in_img_data = self.map_labels(in_img_data,self.labels )


        padding0 = (self.feature_shape[0] + step_size[0] + 1, self.feature_shape[0] + step_size[0] + 1)
        padding1 = (self.feature_shape[1] + step_size[1] + 1, self.feature_shape[1] + step_size[1] + 1)
        padding2 = (self.feature_shape[2] + step_size[2] + 1, self.feature_shape[2] + step_size[2] + 1)




        in_img_data_pad = np.pad(in_img_data, (padding0, padding1, padding2), 'constant', constant_values=0)

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
        if is_label_img == True:
            for patch_iter in range(len(idx_x)):

                curr_patch = in_img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                             idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1],
                             idx_z[patch_iter]:idx_z[patch_iter] + self.feature_shape[2]]

                ## CONVERT THE LABELS HERE.
                # map labels to 0,n_labels

                # curr_patch = self.map_labels(curr_patch, self.labels)
                patches.append(curr_patch)


        else:
            for patch_iter in range(len(idx_x)):
                curr_patch = in_img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                             idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1],
                             idx_z[patch_iter]:idx_z[patch_iter] + self.feature_shape[2]]
                patches.append(curr_patch)


        patches = np.asarray(patches)
        if (is_label_img == False) and  (self.preprocessing == True):
            orig_shape = patches.shape
            patches = patches.reshape(orig_shape[0],-1)
            patches = preprocessing.scale(patches)
            patches = patches.reshape(orig_shape)

        # add channel as a dimension for keras
        newshape = list(patches.shape)
        newshape.append(1)
        # print newshape
        # print newshape.__class__
        patches = np.reshape(patches, newshape)
        indices = np.concatenate((idx_x.reshape(-1,1), idx_y.reshape(-1, 1), idx_z.reshape(-1, 1)), axis=1)



        return patches, np.int32(indices), padded_img_size

    def build_image_from_patches(self, in_patches, indices, padded_img_size, patch_crop_size, step_size):
        ''' patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]'''
        out_img_data = np.zeros(padded_img_size)
        count_img_data = np.zeros(padded_img_size)

        out_feature_shape = in_patches.shape[1:4]



        idx_x = indices[:, 0]
        idx_y = indices[:, 1]
        idx_z = indices[:, 2]

        patch_mask = np.zeros(out_feature_shape)
        patch_mask[0 + patch_crop_size[0]: out_feature_shape[0] - patch_crop_size[0],
        0 + patch_crop_size[1]: out_feature_shape[1] - patch_crop_size[1],
        0 + patch_crop_size[2]: out_feature_shape[2] - patch_crop_size[2]] = 1

        for patch_iter in range(len(idx_x)):
            out_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2]] += \
                np.multiply(np.reshape(in_patches[patch_iter, :],out_feature_shape), patch_mask)

            count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2]] += patch_mask

        out_img_data = np.divide(out_img_data, count_img_data)
        out_img_data[np.isnan(out_img_data)] = 0
        # remove the padding
        unpadded_img_size = padded_img_size - np.multiply(np.asarray(self.feature_shape) + step_size + 1, 2)
        padding = np.asarray(self.feature_shape) + step_size + 1

        print("padding is "+ str(np.multiply(np.asarray(self.feature_shape) + step_size + 1, 2)))
        print("unpadded image size is "+str(unpadded_img_size))

        out_img_data = out_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                       padding[1]:padding[1] + unpadded_img_size[1],
                       padding[2]:padding[2] + unpadded_img_size[2]]
        count_img_data = count_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                         padding[1]:padding[1] + unpadded_img_size[1],
                         padding[2]:padding[2] + unpadded_img_size[2]]

        return out_img_data, count_img_data


    def build_seg_from_patches(self, in_patches, indices, padded_img_size, patch_crop_size, step_size):
        ''' patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]'''
        out_img_data = np.zeros(padded_img_size)
        count_img_data = np.zeros(padded_img_size)

        out_feature_shape = in_patches.shape[1:5]



        idx_x = indices[:, 0]
        idx_y = indices[:, 1]
        idx_z = indices[:, 2]

        patch_mask = np.zeros(out_feature_shape)
        patch_mask[0 + patch_crop_size[0]: out_feature_shape[0] - patch_crop_size[0],
        0 + patch_crop_size[1]: out_feature_shape[1] - patch_crop_size[1],
        0 + patch_crop_size[2]: out_feature_shape[2] - patch_crop_size[2],:] = 1

        for patch_iter in range(len(idx_x)):
            out_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2],:] += \
                np.multiply(np.reshape(in_patches[patch_iter, :],out_feature_shape), patch_mask)

            count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2],:] += patch_mask

        out_img_data = np.divide(out_img_data, count_img_data)
        out_img_data[np.isnan(out_img_data)] = 0
        # remove the padding
        # unpadded_img_size = padded_img_size[0:3] - np.multiply(self.feature_shape, 2)

        unpadded_img_size = padded_img_size[0:3] - np.multiply(np.asarray(self.feature_shape) + step_size + 1, 2)
        padding = np.asarray(self.feature_shape) + step_size + 1

        print("padding is "+ str(np.multiply(np.asarray(self.feature_shape) + step_size + 1, 2)))
        print("unpadded image size is "+str(unpadded_img_size))

        out_img_data = out_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                       padding[1]:padding[1] + unpadded_img_size[1],
                       padding[2]:padding[2] + unpadded_img_size[2]]
        count_img_data = count_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                         padding[1]:padding[1] + unpadded_img_size[1],
                         padding[2]:padding[2] + unpadded_img_size[2]]



        # out_img_data = out_img_data[self.feature_shape[0]:self.feature_shape[0] + unpadded_img_size[0],
        #                self.feature_shape[1]:self.feature_shape[1] + unpadded_img_size[1],
        #                self.feature_shape[2]:self.feature_shape[2] + unpadded_img_size[2], :]
        # count_img_data = count_img_data[self.feature_shape[0]:self.feature_shape[0] + unpadded_img_size[0],
        #                  self.feature_shape[1]:self.feature_shape[1] + unpadded_img_size[1],
        #                  self.feature_shape[2]:self.feature_shape[2] + unpadded_img_size[2], :]

        label_img_data = np.argmax(out_img_data, axis=-1)
        label_img_data = self.map_inv_labels(label_img_data, self.labels)
        return out_img_data, label_img_data, count_img_data


    def generate_src_trg_training_data(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                                       target_label_list = None, step_size=None,
                                       ):
        print('Creating source image patches.')
        if self.use_patches == True:
            (_,_, indices_list) = self.create_training_feature_array(source_filenames, 'src',
                                                                     indices_list=None, step_size=step_size,
                                                                     is_label_img=is_src_label_img)
            print('Creating target image patches.')
            self.create_training_feature_array(target_filenames,  'trg', indices_list,
                                               step_size=step_size,
                                               is_label_img=is_trg_label_img)
        else:
            (_, _, indices_list) = self.create_training_feature_array(source_filenames, 'src',
                                                                      indices_list=None, step_size=None,
                                                                      is_label_img=is_src_label_img)

            self.create_training_label_array(target_label_list ,'trg', indices_list)




    def generate_src_trg_validation_data(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                                       target_label_list=None, step_size=None):
        print('Creating source image patches.')

        if self.use_patches == True:
            (_,_, indices_list) = self.create_training_feature_array(source_filenames,
                                                              'src_validation', indices_list=None,
                                                                     step_size=step_size, is_label_img=is_src_label_img)

            print('Creating target image patches.')
            self.create_training_feature_array(target_filenames,  'trg_validation',
                                              indices_list=indices_list, step_size=step_size,
                                              is_label_img=is_trg_label_img)
        else:
            (_, _, indices_list) = self.create_training_feature_array(source_filenames, 'src_validation',
                                                                      indices_list=None, step_size=None,
                                                                      is_label_img=is_src_label_img)

            self.create_training_label_array(target_label_list ,'trg_validation', indices_list)



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
                y_list.append(self.data_storage.root.trg_validation[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list

    def map_labels(self, input_label_patch, input_label_list):
        output_label_patch = np.zeros(input_label_patch.shape)
        # 0th label is always 0
        for  out_label, in_label in enumerate(input_label_list):
            output_label_patch[input_label_patch == in_label] = out_label

        return output_label_patch

    def map_inv_labels(self, input_label_patch, input_label_list):
        output_label_patch = np.zeros(input_label_patch.shape)
        # 0th label is always 0
        for curr_label in range(len(input_label_list)):
            output_label_patch[input_label_patch == curr_label] = input_label_list[curr_label]

        return output_label_patch



    def get_multi_class_labels(self, data, n_labels):
        """
        Translates a label map into a set of binary labels.
        :param data: numpy array containing the label map with shape: (n_samples, 1, ...).
        :param n_labels: number of labels.
        :param labels: integer values of the labels.
        :return: binary numpy array of shape: (n_samples, n_labels, ...)

        """
        labels = range(n_labels)
        new_shape = list(data.shape[0:4]) + [n_labels] #[data.shape[0], n_labels] + list(data.shape[2:])


        y = np.zeros(new_shape, np.int8)
        for label_index in range(n_labels):
            if labels is not None:
                y[:, :, :, :, label_index][data[:,:,:,:, 0] == labels[label_index]] = 1
            else:
                y[:, :, :, :, label_index][data[:,:,:, :,0] == (label_index + 1)] = 1
        return y

    def convert_data(self, x_list, y_list, n_labels=1, labels=None):
        x = np.asarray(x_list)
        y = np.asarray(y_list)
        if n_labels == 1:
            y[y > 0] = 1
        elif n_labels > 1:
            y = self.get_multi_class_labels(y, n_labels=n_labels)
        return x, y

    def seg_training_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src[index])
                y_list.append(self.data_storage.root.trg[index])


            # probabilistically augment (rotate a batch
            x, y = self.convert_data(x_list, y_list, n_labels=self.n_labels, labels=self.labels)


            yield x, y

    def seg_training_generator_augment(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src[index])
                y_list.append(self.data_storage.root.trg[index])


            # probabilistically augment (rotate a batch
            # x, y = self.convert_data(x_list, y_list, n_labels=self.n_labels, labels=self.labels)
            x = np.asarray(x_list)
            y = np.asarray(y_list)

            yield x, y


    def seg_validation_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src_validation.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src_validation[index])
                y_list.append(self.data_storage.root.trg_validation[index])
            x, y = self.convert_data(x_list, y_list, n_labels=self.n_labels, labels=self.labels)
            yield x, y




    def training_label_generator(self, batch_size):
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

    def validation_label_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src_validation.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src_validation[index])
                y_list.append(self.data_storage.root.trg_validation[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list

if __name__ == "__main__":
    os.environ['LD_LIBRARY_PATH'] = '/usr/pubsw/packages/CUDA/lib64/'

    aseg_labels = np.loadtxt('aseg_labels.txt')
    aparcaseg_labels = np.loadtxt('aparc+aseg_labels.txt')

    feature_shape = (32,32,32)
    f = 32
    d = 2
    curr_unet = DeepImageSynth(unet_num_filters=f, unet_depth=d,
                                              unet_downsampling_factor=1,
                                              feature_shape=feature_shape,
                                              storage_loc="disk",
                                              temp_folder="/local_mount/space/bhim/1/users/aj660/tmp",
                                              n_labels=len(aseg_labels),
                                              labels=list(aseg_labels))


    fs_dir = '/autofs/space/mreuter/users/amod/pfizer_dataset_analysis/data/fs_syn_reg_dir_v1/freesurfer6p0_skullstripped_v1'
    src_img_input_type = 'orig/001'
    src_scanner = 'TRIOmecho'
    src_filenames = fetch_training_data_files(fs_dir, src_scanner, src_img_input_type, np.array([0,2,8,9]))

    src_seg_img_input_type = 'aparc+aseg'
    src_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([0,2,8,9]))

    trg_img_input_type = 'aseg'
    trg_scanner = 'TRIOmprage_1'
    trg_seg_img_input_type = 'aparc+aseg'
    trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([0,2,8,9]))
    trg_seg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_seg_img_input_type, np.array([0,2,8,9]))

    src_img_input_type = 'orig/001'
    src_scanner = 'TRIOmecho'
    val_src_filenames = fetch_training_data_files(fs_dir, src_scanner, src_img_input_type, np.array([3]))
    src_seg_img_input_type = 'aparc+aseg'
    val_src_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([3]))

    trg_img_input_type = 'aseg'
    trg_scanner = 'TRIOmprage_1'
    trg_seg_img_input_type = 'aparc+aseg'
    val_trg_filenames = fetch_training_data_files(fs_dir, trg_scanner, trg_img_input_type, np.array([3]))
    val_trg_seg_filenames = fetch_training_data_files(fs_dir, src_scanner, src_seg_img_input_type, np.array([3]))


    step_size = (4,4,4)

    curr_unet.load_training_images(src_filenames, trg_filenames,
                                   src_seg_filenames, trg_seg_filenames, step_size=step_size)

    curr_unet.load_validation_images(val_src_filenames, val_trg_filenames,
                                     val_src_seg_filenames, val_trg_seg_filenames, step_size=step_size)

    # output_dir = "/autofs/space/mreuter/users/amod/deep_learn/results/seg_unet_depth_"+ \
    #              str(d) +"_filts_" + str(f) +"_patch"+str(feature_shape[0])+"x"+str(feature_shape[1])+ \
    #              "x"+str(feature_shape[2])
    # subprocess.call(['mkdir', '-p', output_dir])
    #
    # test_src_filenames = fetch_training_data_files(fs_dir, src_scanner, src_img_input_type, np.array([1,4,5,6,7,10,11,12]))
    # out_dir = opj(output_dir, src_scanner)
    # subprocess.call(['mkdir', '-p', out_dir])
    # print(curr_unet.model.summary())
    # curr_unet.train_network(output_prefix=opj(output_dir, src_scanner+"_to_"+trg_scanner), epochs=5,
    #                         initial_epoch=1,batch_size=32, steps_per_epoch=10000,
    #                         save_per_epoch=True, save_weights=True)
    #

    # fg.generate_src_trg_validation_data(val_src_filenames, val_trg_filenames, val_src_seg_filenames, val_trg_seg_filenames, step_size=[16,16,16])
