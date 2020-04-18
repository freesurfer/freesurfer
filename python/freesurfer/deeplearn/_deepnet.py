import os
import sys
import numpy as np
import nibabel as nib
import glob
import tables
import pickle
import time
import tempfile
import scipy
import random
import sklearn
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import ReduceLROnPlateau, TensorBoard, ModelCheckpoint
from tensorflow.keras import backend
from .. import wm_peak_normalize, robust_normalize, wm_peak_normalize_t2w
from . import *


class MultiGPUCheckpointCallback(tf.keras.callbacks.Callback):
    def __init__(self, output_prefix, save_per_epoch, save_weights, initial_epoch=1):
        super(MultiGPUCheckpointCallback, self).__init__()
        self.output_prefix = output_prefix
        self.save_per_epoch = save_per_epoch
        self.save_weights = save_weights
        self.initial_epoch = initial_epoch

    def detach_model(m):
        for l in m.layers:
            if l.name == 'model_1': return l
        return m

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


class MRIDeepNetCallback(tf.keras.callbacks.Callback):
    def __init__(self, output_prefix, save_per_epoch, save_weights, initial_epoch=1):
        self.output_prefix = output_prefix
        self.save_per_epoch = save_per_epoch
        self.save_weights = save_weights
        self.initial_epoch = initial_epoch
        super(MRIDeepNetCallback, self).__init__()

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


class MRIDeepNet(object):
    def __init__(self, unet_num_filters, unet_depth, unet_downsampling_factor, feature_shape, depth_per_level=2, dim=3,
                 kernel_size=(3, 3, 3),
                 storage_loc="memory", temp_folder=os.getcwd(), num_input_channels=1, channel_names=['t1w'],
                 out_channel_names=['seg'],
                 n_labels=0, labels=[], net='unet', loss='mean_absolute_error',
                 initial_learning_rate=0.00001, use_patches=True, wmp_standardize=True, rob_standardize=True,
                 fcn=True, num_gpus=1, preprocessing=False, augment=False, num_outputs=1,
                 nmr_augment=False, use_tal=True, rare_label_list=[], use_slices=False, orientation=None,
                 add_modality_channel=False):
        self.net = net
        self.num_input_channels = num_input_channels
        self.channel_names = channel_names
        self.out_channel_names = out_channel_names

        self.feature_shape = feature_shape  # should already include the num_input_channels (32,32,32,4)
        self.storage_loc = storage_loc
        self.temp_folder = temp_folder
        self.wmp_standardize = wmp_standardize
        self.rob_standardize = rob_standardize,
        self.fcn = fcn  # fully convolutional or not boolean
        self.num_gpus = num_gpus
        self.preprocessing = preprocessing
        self.augment = augment
        self.nmr_augment = nmr_augment
        self.use_tal = use_tal
        self.rare_label_list = rare_label_list

        self.dim = len(feature_shape) - 1
        self.num_outputs = num_outputs

        self.orientation = orientation
        self.kernel_size = kernel_size
        self.depth_per_level = depth_per_level

        self.use_slices = use_slices
        self.add_modality_channel = add_modality_channel

        if net == 'unet':
            self.unet_downsampling_factor = unet_downsampling_factor
            self.unet_num_filters = unet_num_filters
            self.unet_depth = unet_depth
            self.model, self.parallel_model = unet_model_3d(feature_shape, num_filters=unet_num_filters,
                                                            unet_depth=unet_depth,
                                                            downsize_filters_factor=unet_downsampling_factor,
                                                            pool_size=(2, 2, 2),
                                                            n_labels=n_labels, loss=loss,
                                                            initial_learning_rate=initial_learning_rate,
                                                            deconvolution=False, use_patches=use_patches,
                                                            num_gpus=num_gpus)
        elif net == 'resnet':
            self.unet_downsampling_factor = unet_downsampling_factor
            self.unet_num_filters = unet_num_filters
            self.unet_depth = unet_depth
            self.model, self.parallel_model = unet_model(feature_shape, num_filters=unet_num_filters,
                                                         unet_depth=unet_depth,
                                                         downsize_filters_factor=unet_downsampling_factor,
                                                         pool_size=(2, 2, 2),
                                                         n_labels=n_labels, loss=loss,
                                                         initial_learning_rate=initial_learning_rate,
                                                         deconvolution=False, use_patches=use_patches,
                                                         num_gpus=num_gpus)


        elif net == 'atrousnet':
            self.model = atrous_net(feature_shape, unet_num_filters, initial_learning_rate=initial_learning_rate,
                                    loss='mean_absolute_error')
        elif net == 'noisenet':
            self.model = noise_net(feature_shape, unet_num_filters, initial_learning_rate=initial_learning_rate,
                                   loss='mean_absolute_error')
        elif (net == 'class_net') & (fcn == False):
            # i.e. not a FCN like the above. A single label for the input feature, not a patch of labels as above
            self.model = class_net(feature_shape=feature_shape, dim=dim, unet_num_filters=unet_num_filters,
                                   n_labels=n_labels, initial_learning_rate=initial_learning_rate, loss=loss)
        elif net == 'unet2d':
            self.model, self.parallel_model = unet_model_2d(feature_shape, num_filters=unet_num_filters,
                                                            unet_depth=unet_depth,
                                                            downsize_filters_factor=unet_downsampling_factor,
                                                            pool_size=(2, 2),
                                                            n_labels=n_labels, loss=loss,
                                                            initial_learning_rate=initial_learning_rate,
                                                            deconvolution=False, use_patches=use_patches,
                                                            num_gpus=num_gpus, num_outputs=num_outputs)

        self.model_trained = False
        self.model_compiled = False
        self.labels = labels
        self.n_labels = len(labels)
        self.feature_generator = FeatureGenerator(self.feature_shape, temp_folder=temp_folder,
                                                  storage_loc=storage_loc, labels=labels, n_labels=n_labels,
                                                  wmp_standardize=wmp_standardize, rob_standardize=rob_standardize,
                                                  use_patches=use_patches, dim=dim,
                                                  preprocessing=self.preprocessing, augment=augment,
                                                  nmr_augment=nmr_augment,
                                                  num_outputs=num_outputs,
                                                  num_input_channels=num_input_channels, channel_names=channel_names,
                                                  out_channel_names=out_channel_names,
                                                  use_tal=use_tal,
                                                  rare_label_list=rare_label_list,
                                                  orientation=orientation,
                                                  use_slices=use_slices,

                                                  )
        self._weight_file = None

    @classmethod
    def from_file(cls, model_filename, loss, storage_loc='memory', temp_folder=os.getcwd(),
                  net='unet', n_labels=0, labels=[], num_input_channels=1, channel_names=['t1w'],
                  out_channel_names=['seg'], initial_learning_rate=0.00001, use_patches=True,
                  wmp_standardize=True, rob_standardize=True, fcn=True, num_gpus=1,
                  preprocessing=False, augment=False, num_outputs=1,
                  nmr_augment=False, use_tal=True, rare_label_list=[], use_slices=False,
                  add_modality_channel=False,

                  ):
        if loss == 'dice_coef_loss':
            model = load_model(model_filename, custom_objects={'dice_coef_loss': dice_coef_loss})
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

        kernel_size = layer_list[1]['config']['kernel_size']

        depth_per_level = 0
        for iter in range(len(layer_list)):
            layer_name = layer_list[iter]['class_name']
            if 'Conv' in layer_name:
                depth_per_level = depth_per_level + 1
            elif 'Pool' in layer_name:
                break

        num_pools = 0
        for layer in layer_list:
            if ((layer['class_name'] == 'MaxPooling3D') | ((layer['class_name'] == 'MaxPooling2D'))):
                num_pools = num_pools + 1

        unet_depth = num_pools + 1  # 2 * num_pools - 1
        print(unet_depth)

        feature_shape = tuple(input_shape[1:])

        if net == 'unet':

            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                           unet_downsampling_factor=unet_downsampling_factor,
                           feature_shape=feature_shape, loss=loss, storage_loc=storage_loc,
                           n_labels=n_labels, labels=labels,
                           temp_folder=temp_folder, net='unet',
                           num_input_channels=num_input_channels, channel_names=channel_names,
                           out_channel_names=out_channel_names, initial_learning_rate=initial_learning_rate,
                           use_patches=use_patches,
                           wmp_standardize=wmp_standardize, rob_standardize=rob_standardize, fcn=fcn, num_gpus=1,
                           preprocessing=preprocessing, augment=augment, num_outputs=num_outputs,
                           nmr_augment=nmr_augment, use_tal=use_tal, rare_label_list=[]

                           )
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True
            print('Loaded model file: %s' % model_filename)
        elif net == 'atrousnet':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth,
                           unet_downsampling_factor=unet_downsampling_factor,
                           feature_shape=feature_shape, storage_loc=storage_loc, temp_folder=temp_folder,
                           net='atrousnet')
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
        elif net == 'unet_2d_v1':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth, depth_per_level=depth_per_level,
                           unet_downsampling_factor=unet_downsampling_factor, kernel_size=kernel_size,
                           feature_shape=feature_shape, loss=loss, storage_loc=storage_loc,
                           n_labels=n_labels, labels=labels,
                           temp_folder=temp_folder, net='unet_2d_v1',
                           wmp_standardize=wmp_standardize, rob_standardize=rob_standardize,
                           use_tal=use_tal, use_slices=use_slices, use_patches=use_patches,
                           add_modality_channel=add_modality_channel)
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True
            print('Loaded model file: %s' % model_filename)

        return cls_init

    @classmethod
    def from_file_old_model(cls, model_filename, loss, storage_loc='memory', temp_folder=os.getcwd(),
                            net='class_net', n_labels=0):
        model = tf.keras.load_model(model_filename)

        input_shape = model.input_shape
        unet_num_filters = model.get_config()[0]['config']['filters']
        dim = len(input_shape) - 2

        unet_depth = len(model.get_config()[0]) / 4

        feature_shape = tuple(input_shape[1:-1])
        if net == 'class_net':
            cls_init = cls(unet_num_filters=unet_num_filters, unet_depth=unet_depth, unet_downsampling_factor=1,
                           feature_shape=feature_shape, dim=dim, storage_loc=storage_loc,
                           temp_folder=temp_folder, n_labels=n_labels,
                           net='class_net')
            cls_init.model = model
            cls_init.model_trained = True
            cls_init.model_compiled = True

        return cls_init

    def load_training_images(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                             source_seg_filenames=None, target_seg_filenames=None, step_size=None, preprocessing=False):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_data(source_filenames, target_filenames,
                                                              is_src_label_img, is_trg_label_img,
                                                              step_size=step_size, )

    def load_training_image_slices(self, source_filenames, target_filenames, is_src_label_img,
                                   is_trg_label_img, orientation, ):
        if self.model is None:
            raise RuntimeError('Model does not exist')
        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_slice_data(source_filenames, target_filenames,
                                                                    is_src_label_img, is_trg_label_img, orientation)

    def load_validation_image_slices(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                                     orientation, ):
        if self.model is None:
            raise RuntimeError('Model does not exist')
        self.feature_generator.generate_src_trg_validation_slice_data(source_filenames, target_filenames,
                                                                      is_src_label_img, is_trg_label_img, orientation)

    def load_training_slices_and_labels(self, source_filenames, label_list, is_src_label_img):
        if self.model is None:
            raise RuntimeError('Model does not exist')
        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_training_data(source_filenames, target_filenames=None,
                                                              is_trg_label_img=False,
                                                              is_src_label_img=is_src_label_img,
                                                              target_label_list=label_list)

    def load_validation_images(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                               source_seg_filenames=None, target_seg_filenames=None, step_size=None,
                               preprocessing=False):
        if self.model is None:
            raise RuntimeError('Model does not exist')
        self.feature_generator.generate_src_trg_validation_data(source_filenames, target_filenames,
                                                                is_src_label_img, is_trg_label_img,
                                                                step_size=step_size)

    def load_validation_slices_and_labels(self, source_filenames, label_list, is_src_label_img):
        if self.model is None:
            raise RuntimeError('Model does not exist')

        self.feature_generator.generate_src_trg_validation_data(source_filenames, target_filenames=None,
                                                                is_trg_label_img=False,
                                                                is_src_label_img=is_src_label_img,
                                                                target_label_list=label_list)

    def load_training_validation_image_slices(self, source_filename, target_filename, step_size=None):
        # for the oct mgz file
        self.feature_generator.create_data_storage()
        self.feature_generator.generate_src_trg_slice_training_validation(source_filename, target_filename,
                                                                          step_size=step_size)

    def train_network(self, output_prefix, epochs=5, initial_epoch=1, batch_size=64, steps_per_epoch=10000,
                      validation_steps=1000,
                      optimizer='adam', save_per_epoch=False, save_weights=True, finetune=False, focus='ALL'):
        print('Beginning Training. Using %s backend with "%s" data format.' % (backend.backend(),
                                                                               backend.image_data_format()))
        if self._weight_file is not None:
            self.model.load_weights(self._weight_file)

        reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5, min_lr=0.000001)
        # modelcp = ModelCheckpoint(output_prefix, monitor='val_loss', verbose=0, save_best_only=False,
        #                                 save_weights_only=False, mode='auto', period=1)

        print('Training model...')
        if self.n_labels == 0:
            if self.num_gpus == 1:
                callback = MRIDeepNetCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)
                if len(self.out_channel_names) == 1:
                    out_channel_name = self.out_channel_names[0]
                    if out_channel_name == 't1beta':
                        self.model.fit_generator(
                            generator=self.feature_generator.training_generator_t1beta(
                                batch_size=batch_size * self.num_gpus),
                            epochs=epochs,
                            validation_data=self.feature_generator.validation_generator_t1beta(
                                batch_size=batch_size * self.num_gpus),
                            validation_steps=validation_steps,
                            steps_per_epoch=steps_per_epoch,
                            callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)
                    elif out_channel_name == 'pdbeta':
                        self.model.fit_generator(
                            generator=self.feature_generator.training_generator_pdbeta(
                                batch_size=batch_size * self.num_gpus),
                            epochs=epochs,
                            validation_data=self.feature_generator.validation_generator_pdbeta(
                                batch_size=batch_size * self.num_gpus),
                            validation_steps=validation_steps,
                            steps_per_epoch=steps_per_epoch,
                            callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)

                    elif out_channel_name == 't2beta':
                        self.model.fit_generator(
                            generator=self.feature_generator.training_generator_t2beta(
                                batch_size=batch_size * self.num_gpus),
                            epochs=epochs,
                            validation_data=self.feature_generator.validation_generator_t2beta(
                                batch_size=batch_size * self.num_gpus),
                            validation_steps=validation_steps,
                            steps_per_epoch=steps_per_epoch,
                            callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)

                    else:
                        self.model.fit_generator(
                            generator=self.feature_generator.training_generator(batch_size=batch_size * self.num_gpus),
                            epochs=epochs,
                            validation_data=self.feature_generator.validation_generator(
                                batch_size=batch_size * self.num_gpus),
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

        elif (self.n_labels > 1) & (self.fcn == True):
            if self.num_gpus == 1:
                callback = MRIDeepNetCallback(output_prefix, save_per_epoch, save_weights, initial_epoch)
                if ((self.nmr_augment == True) & (self.use_tal == True)):
                    # self.model.fit_generator(generator=self.feature_generator.seg_training_generator_multichannel_nmr_t2(batch_size=batch_size*self.num_gpus),
                    #                  epochs=epochs,
                    #                  validation_data=self.feature_generator.seg_validation_generator_multichannel_nmr(batch_size=batch_size*self.num_gpus),
                    #                  validation_steps=validation_steps,
                    #                  steps_per_epoch=steps_per_epoch,
                    #                  callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=10)

                    self.model.fit_generator(
                        generator=self.feature_generator.dynamic_seg_training_generator_multichannel_nmr_t2(
                            batch_size=batch_size * self.num_gpus,
                            focus=focus),
                        epochs=epochs,
                        validation_data=self.feature_generator.dynamic_seg_validation_generator_multichannel_nmr_t2(
                            batch_size=batch_size * self.num_gpus),
                        validation_steps=validation_steps,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=10)
                elif ((self.nmr_augment == True) & (self.use_tal == True) & (finetune == True)):
                    self.model.fit_generator(
                        generator=self.feature_generator.seg_training_generator_multichannel_nmr_finetune(
                            batch_size=batch_size * self.num_gpus),
                        epochs=epochs,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=10)

                elif ((self.nmr_augment == True) & (self.use_tal == False) & (self.use_slices == False)):
                    self.model.fit_generator(
                        generator=self.feature_generator.dynamic_seg_training_generator_singlechannel_nmr(
                            batch_size=batch_size * self.num_gpus, focus=focus),
                        epochs=epochs,
                        validation_data=self.feature_generator.dynamic_seg_validation_generator_singlechannel_nmr(
                            batch_size=batch_size * self.num_gpus),
                        validation_steps=validation_steps,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)


                elif ((self.nmr_augment == True) & (self.use_tal == False) & (self.use_slices == True)):
                    print('focus is ' + focus)
                    self.model.fit_generator(
                        generator=self.feature_generator.dynamic_seg_training_generator_singlechannel_nmr_slice(
                            batch_size=batch_size * self.num_gpus, focus=focus),
                        epochs=epochs,
                        validation_data=self.feature_generator.dynamic_seg_validation_generator_singlechannel_nmr_slice(
                            batch_size=batch_size * self.num_gpus),
                        validation_steps=validation_steps,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)



                elif ((self.nmr_augment == False) & (self.use_tal == False)):
                    self.model.fit_generator(
                        generator=self.feature_generator.dynamic_seg_training_generator_singlechannel(
                            batch_size=batch_size * self.num_gpus),
                        epochs=epochs,
                        validation_data=self.feature_generator.dynamic_seg_validation_generator_singlechannel(
                            batch_size=batch_size * self.num_gpus),
                        validation_steps=validation_steps,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)
                else:
                    self.model.fit_generator(
                        generator=self.feature_generator.seg_training_generator_multichannel(
                            batch_size=batch_size * self.num_gpus),
                        epochs=epochs,
                        validation_data=self.feature_generator.seg_validation_generator_multichannel(
                            batch_size=batch_size * self.num_gpus),
                        validation_steps=validation_steps,
                        steps_per_epoch=steps_per_epoch,
                        callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)
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
            self.model.fit_generator(
                generator=self.feature_generator.training_label_generator(batch_size=batch_size * self.num_gpus),
                epochs=epochs,
                validation_data=self.feature_generator.validation_label_generator(
                    batch_size=batch_size * self.num_gpus),
                validation_steps=validation_steps,
                steps_per_epoch=steps_per_epoch,
                callbacks=[callback, reduce_lr, ], verbose=1, max_queue_size=100)

    def synthesize_image(self, in_img_file_list, channel_names, out_img_filename, step_size, scale_output=1):
        num_channels = len(in_img_file_list)
        if num_channels != len(channel_names):
            print('Channel number and input number of files do not match!')
        input_features = list()
        for iter_channel in range(num_channels):
            in_img_file = in_img_file_list[iter_channel]
            in_img = nib.load(in_img_file)
            if iter_channel == 0:
                (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file,
                                                                                                   intensity_threshold=0,
                                                                                                   step_size=step_size,
                                                                                                   is_label_img=False,
                                                                                                   indices=None,
                                                                                                   channel_name=
                                                                                                   channel_names[
                                                                                                       iter_channel])
                input_features.append(in_patches)
            else:
                (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file,
                                                                                                   intensity_threshold=0,
                                                                                                   step_size=step_size,
                                                                                                   is_label_img=False,
                                                                                                   indices=in_indices,
                                                                                                   channel_name=
                                                                                                   channel_names[
                                                                                                       iter_channel])
                input_features.append(in_patches)

        input_feature_array = np.concatenate(input_features, axis=4)
        print(input_feature_array.shape)
        out_patches = self.model.predict(input_feature_array)

        # in_img = nib.load(in_img_file)
        # print("Shape is " + str(in_img.get_data().shape))
        # (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file, intensity_threshold=0,
        #                                                                  step_size=step_size, is_label_img=False, indices=None)

        # if self.num_gpus > 1:
        #     out_patches = self.parallel_model.predict(in_patches)
        # else:
        #     out_patches = self.model.predict(in_patches)

        patch_crop_size = [1, 1, 1]  # should be filter_size/2

        print("padded image size " + str(padded_img_size))
        out_img_data, count_img_data = self.feature_generator.build_image_from_patches(out_patches, in_indices,
                                                                                       padded_img_size, patch_crop_size,
                                                                                       step_size)

        out_img_data = out_img_data * scale_output
        in_img.set_data_dtype(np.dtype('float32'))

        print("Out data shape is " + str(out_img_data.shape))

        out_img = nib.MGHImage(out_img_data, in_img.affine, in_img.header)
        nib.save(out_img, out_img_filename)

    def predict_slice_segmentation(self, in_img_file_list, channel_names, orientation, out_soft_filename,
                                   out_hard_filename):
        num_channels = len(in_img_file_list)
        # these should match the channel_names lenght
        if num_channels != len(channel_names):
            print('Channel number and input number of files do not match!')
        input_features = list()
        for iter_channel in range(num_channels):
            in_img_file = in_img_file_list[iter_channel]
            in_img = nib.load(in_img_file)
            if iter_channel == 0:
                (in_slices, in_indices) = self.feature_generator.extract_slices(in_img_file, intensity_threshold=0,
                                                                                indices=None,
                                                                                is_label_img=False,
                                                                                orientation=orientation)
                input_features.append(in_slices)
            else:
                (in_slices, in_indices) = self.feature_generator.extract_slices(in_img_file, intensity_threshold=0,
                                                                                indices=in_indices,
                                                                                is_label_img=False,
                                                                                orientation=orientation)
                input_features.append(in_slices)

        input_feature_array = np.concatenate(input_features, axis=-1)
        print(input_feature_array.shape)

        if self.num_gpus > 1:
            out_slices = self.parallel_model.predict(input_feature_array)
        else:
            out_slices = self.model.predict(input_feature_array)

        label_slices = np.argmax(out_slices, axis=-1)
        label_slices = self.feature_generator.map_inv_labels(label_slices, self.labels)

        out_slices = np.moveaxis(out_slices, 0, -1)
        label_slices = np.moveaxis(label_slices, 0, -1)

        out_img = nib.MGHImage(out_slices, in_img.affine, in_img.header)
        nib.save(out_img, out_soft_filename)

        label_img = nib.MGHImage(label_slices, in_img.affine, in_img.header)
        nib.save(label_img, out_hard_filename)

    def predict_segmentation(self, in_img_file_list, channel_names, out_soft_filename, out_hard_filename, step_size):
        num_channels = len(in_img_file_list)
        # these should match the channel_names lenght
        if num_channels != len(channel_names):
            print('Channel number and input number of files do not match!')
        input_features = list()
        for iter_channel in range(num_channels):
            in_img_file = in_img_file_list[iter_channel]
            in_img = nib.load(in_img_file)
            if iter_channel == 0:
                (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file,
                                                                                                   intensity_threshold=0,
                                                                                                   step_size=step_size,
                                                                                                   is_label_img=False,
                                                                                                   indices=None,
                                                                                                   channel_name=
                                                                                                   channel_names[
                                                                                                       iter_channel])
                input_features.append(in_patches)
            else:
                (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file,
                                                                                                   intensity_threshold=0,
                                                                                                   step_size=step_size,
                                                                                                   is_label_img=False,
                                                                                                   indices=in_indices,
                                                                                                   channel_name=
                                                                                                   channel_names[
                                                                                                       iter_channel])
                input_features.append(in_patches)

        input_feature_array = np.concatenate(input_features, axis=4)
        # input_feature_array = input_feature_array.reshape(input_feature_array.shape[0:-1])
        print(input_feature_array.shape)

        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(input_feature_array)
        else:
            out_patches = self.model.predict(input_feature_array)

        num_labels = out_patches.shape[-1]
        padded_img_size_multiple_labels = padded_img_size + (num_labels,)

        patch_crop_size = [1, 1, 1]  # should be filter_size/2

        out_img_data, label_img_data, count_img_data = self.feature_generator.build_seg_from_patches(out_patches,
                                                                                                     in_indices,
                                                                                                     padded_img_size_multiple_labels,
                                                                                                     patch_crop_size,
                                                                                                     step_size)

        out_img = nib.MGHImage(out_img_data, in_img.affine, in_img.header)
        nib.save(out_img, out_soft_filename)

        label_img = nib.MGHImage(label_img_data, in_img.affine, in_img.header)
        nib.save(label_img, out_hard_filename)

    def predict_labels(self, in_img_file):
        in_img = nib.load(in_img_file)
        (in_slices, in_indices) = self.feature_generator.extract_slices(in_img_file, intensity_threshold=0,
                                                                        is_label_img=False, indices=None)
        out_prob = self.model.predict(in_slices)
        out_labels = np.zeros(out_prob.shape)
        out_labels[out_prob > 0.5] = 1
        return out_prob, out_labels

    def predict_slices(self, in_img_file, truth_file):
        imgs = nib.load(in_img_file)
        vecs = nib.load(truth_file)

        img_data = imgs.get_data()
        img_data = img_data[:, :, 0]

        vec_data = vecs.get_data()
        vec_data = vec_data[:, :, 0, :]

        padding0 = (self.feature_shape[0] + step_size[0] + 1, self.feature_shape[0] + step_size[0] + 1)
        padding1 = (self.feature_shape[1] + step_size[1] + 1, self.feature_shape[1] + step_size[1] + 1)
        padding2 = (0, 0)
        img_data_pad = np.pad(img_data, (padding0, padding1), 'constant', constant_values=0)
        vec_data_pad = np.pad(vec_data, (padding0, padding1, padding2), 'constant', constant_values=0)

        padded_img_size = img_data_pad.shape
        intensity_threshold = 0

        (idx_x_fg, idx_y_fg, idx_z_fg) = np.where(vec_data_pad > intensity_threshold)
        min_idx_x_fg = np.min(idx_x_fg) - step_size[0]
        max_idx_x_fg = np.max(idx_x_fg) + step_size[0]
        min_idx_y_fg = np.min(idx_y_fg) - step_size[1]
        max_idx_y_fg = np.max(idx_y_fg) + step_size[1]

        sampled_x = np.arange(min_idx_x_fg, max_idx_x_fg, step_size[0])
        sampled_y = np.arange(min_idx_y_fg, max_idx_y_fg, step_size[1])

        idx_x, idx_y = np.meshgrid(sampled_x, sampled_y, sparse=False, indexing='ij')
        idx_x = idx_x.flatten()
        idx_y = idx_y.flatten()
        # idx_z = idx_z.flatten()

        patches = []
        trg_patches = []
        indices_x = []
        indices_y = []
        # take half the patches (because we have no memory)
        # print('dataset is for ' + dataset)
        # if dataset == 'train':
        #
        # elif dataset == 'val':
        #     iterrange = range((len(idx_x)//3 + 1), (len(idx_x) // 3) + (len(idx_x)//5))

        # print('iter range is ' + str(iterrange[0]) + '_' + str(iterrange[-1]))
        iterrange = range(len(idx_x) // 2, len(idx_x))
        for patch_iter in range(idx_x):
            curr_patch = img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                         idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1]]
            vec_patch = vec_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                        idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1], 0:2]

            if vec_patch.mean() != 0:
                print(patch_iter * 100.0 / len(idx_x))
                patches.append(curr_patch)
                trg_patches.append(vec_patch)
                indices_x.append(idx_x[patch_iter])
                indices_y.append(idx_y[patch_iter])
                # save the patches in out.mgz files
                # currshape = list(curr_patch.shape)
                # currshape.append(1)
                # vecshape = list(vec_patch.shape[0:2])
                # vecshape.append(1)
                # vecshape.append(2)

        patches = np.asarray(patches)
        newshape = list(patches.shape)
        newshape.append(1)
        out_vec = self.model.predict(patches, batch_size=32)
        return patches, out_vec, trg_patches

    def apply_encoder(self, in_img_file, layer_name, step_size):
        encoder_model = tf.keras.Model(inputs=self.model.input, outputs=self.model.get_layer(layer_name).output)
        in_img = nib.load(in_img_file)
        print("Shape is " + str(in_img.get_data().shape))
        (in_patches, in_indices, padded_img_size) = self.feature_generator.extract_patches(in_img_file,
                                                                                           intensity_threshold=0,
                                                                                           step_size=step_size,
                                                                                           is_label_img=False,
                                                                                           indices=None)
        if self.num_gpus > 1:
            out_patches = self.parallel_model.predict(in_patches)
        else:
            out_patches = encoder_model.predict(in_patches)
        return out_patches, in_patches, in_indices, padded_img_size

    def apply_decoder(self, in_patches, layer_name, step_size):
        decoder_model = tf.keras.Model(inputs=self.model.get_layer(layer_name), outputs=self.model.output)
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
                 wmp_standardize=True, rob_standardize=True, use_patches=False, dim=3, preprocessing=False,
                 augment=True,
                 nmr_augment=False,
                 prob_rot_augmentation=0.25, prob_int_augmentation=0.5, num_outputs=1,
                 num_input_channels=1, channel_names=['t1w'], out_channel_names=['trg'],
                 input_normalization_key=['wmp'], output_normalization_key=None, use_tal=True, rare_label_list=None,
                 use_slices=False,
                 orientation=None):

        self.feature_shape = feature_shape

        self.temp_folder = temp_folder
        self.storage_loc = storage_loc

        self.n_labels = n_labels
        self.labels = labels
        self.num_outputs = num_outputs

        self.wmp_standardize = wmp_standardize
        self.rob_standardize = rob_standardize

        self.use_patches = use_patches
        self.dim = dim
        self.preprocessing = preprocessing
        self.augment = augment
        self.nmr_augment = nmr_augment

        self.num_input_channels = num_input_channels
        self.channel_names = channel_names
        self.out_channel_names = out_channel_names
        self.use_tal = use_tal
        self.rare_label_list = rare_label_list
        self.use_slices = use_slices
        self.orientation = orientation

        # dicts to store images in memory
        self.trg_image_dict = {}
        self.trg_fg_indices_dict = {}
        self.trg_fg_labels_dict = {}
        self.trg_slice_indices_dict = {}

        self.val_trg_image_dict = {}
        self.val_trg_fg_indices_dict = {}
        self.val_trg_fg_labels_dict = {}
        self.val_trg_slice_indices_dict = {}

        self.src_image_dict = {}
        self.val_src_image_dict = {}

        for ch in self.out_channel_names:
            self.trg_image_dict[ch] = list()
            self.val_trg_image_dict[ch] = list()

            if self.use_slices == True:
                self.trg_slice_indices_dict[ch] = list()
                self.val_trg_slice_indices_dict[ch] = list()

            else:
                self.trg_fg_indices_dict[ch] = list()
                self.trg_fg_labels_dict[ch] = list()
                self.val_trg_fg_indices_dict[ch] = list()
                self.val_trg_fg_labels_dict[ch] = list()

        for ch in self.channel_names:
            self.src_image_dict[ch] = list()
            self.val_src_image_dict[ch] = list()

        if self.augment == True:
            self.prob_rot_augmentation = prob_rot_augmentation
            self.prob_int_augmentation = prob_int_augmentation
            augment_cubic_transforms_dict = pickle.load(open('pfizer_data_cubic_fits.pkl', 'rb'))
            self.augment_cubic_transforms = np.vstack((augment_cubic_transforms_dict['triomprage'],
                                                       augment_cubic_transforms_dict['triomecho'],
                                                       augment_cubic_transforms_dict['gemprage'],
                                                       augment_cubic_transforms_dict['sonatamecho']))

        self.data_storage = None
        self.storage_filename = None

    def create_data_storage(self):
        if self.storage_loc == 'memory':
            self.data_storage = tables.open_file('tmp_data_storage.h5', 'w', driver='H5FD_CORE',
                                                 driver_core_backing_store=False)
        elif self.storage_loc == 'disk':
            tmp_fp = tempfile.NamedTemporaryFile('w', suffix='.h5', dir=self.temp_folder, delete=False)
            self.storage_filename = tmp_fp.name
            tmp_fp.close()
            self.data_storage = tables.open_file(self.storage_filename, 'w')
        else:
            raise RuntimeError('Choose one of {memory, disk} for storage_loc')

    def load_data_storage(self, filename):
        if self.storage_loc == 'memory':
            self.data_storage = tables.open_file(filename, 'r', driver='core', driver_core_backing_store=False)
        elif self.storage_loc == 'disk':
            self.data_storage = tables.open_file(filename, 'r')
        else:
            raise RuntimeError('Choose one of {memory, disk} for storage_loc')

    def close_data_storage(self):
        print('Closing file')
        self.data_storage.close()

    def generate_src_trg_training_collection(self, source_filenames, target_filenames, is_src_label_img,
                                             is_trg_label_img,
                                             target_label_list=None, step_size=None, rare_label_list=None,
                                             ):
        num_subjects = len(target_filenames[0])
        num_source_channels = self.num_input_channels
        num_target_channels = len(self.out_channel_names)

        for ch_iter in range(num_target_channels):
            for subj_iter in range(num_subjects):

                trg_img = nib.load(target_filenames[ch_iter][subj_iter])
                trg_img_data = trg_img.get_data()
                trg_img_data = self.preprocess_image(trg_img_data, self.out_channel_names[ch_iter])

                if self.use_slices == True:
                    if self.orientation is None:
                        print('orientation of slices is not specified')
                        return
                    elif self.orientation == 'coronal':
                        slice_idxs = list()
                        for sl_idx in range(0, trg_img_data.shape[2]):
                            slice = trg_img_data[:, :, sl_idx]
                            nz_idxs = np.argwhere(slice > 0)
                            if nz_idxs.shape[0] > 0:
                                slice_idxs.append(sl_idx)

                    slice_idxs = np.asarray(slice_idxs)
                    trg_img_data = self.map_labels(trg_img_data, self.labels)
                    curr_channel = self.out_channel_names[ch_iter]
                    self.trg_image_dict[curr_channel].append(trg_img_data)
                    self.trg_slice_indices_dict[curr_channel].append(slice_idxs)
                    print(trg_img_data.shape)
                    print(slice_idxs.shape)


                else:

                    if rare_label_list is None:
                        fg_idxs = np.argwhere(trg_img_data > 0)
                        trg_fg = trg_img_data[trg_img_data > 0].flatten()
                        trg_fg = np.reshape(trg_fg, (trg_fg.shape[0], 1))
                    else:
                        print('rare label list is ' + str(rare_label_list))
                        mask = np.isin(trg_img_data, rare_label_list)
                        fg_idxs = np.argwhere(mask == True)
                        trg_fg = trg_img_data[mask == True].flatten()
                        trg_fg = np.reshape(trg_fg, (trg_fg.shape[0], 1))

                    trg_img_data = self.map_labels(trg_img_data, self.labels)
                    print(trg_img_data.shape)
                    print(fg_idxs.shape)
                    curr_channel = self.out_channel_names[ch_iter]
                    self.trg_image_dict[curr_channel].append(trg_img_data)
                    self.trg_fg_indices_dict[curr_channel].append(fg_idxs)
                    self.trg_fg_labels_dict[curr_channel].append(trg_fg)

        for ch_iter in range(num_source_channels):

            curr_channel = self.channel_names[ch_iter]

            for subj_iter in range(num_subjects):
                src_img = nib.load(source_filenames[ch_iter][subj_iter])
                src_img_data = src_img.get_data()
                src_img_data = self.preprocess_image(src_img_data, channel_name=curr_channel)
                self.src_image_dict[curr_channel].append(src_img_data)

    def generate_src_trg_validation_collection(self, source_filenames, target_filenames, is_src_label_img,
                                               is_trg_label_img,
                                               target_label_list=None, step_size=None,
                                               ):
        num_subjects = len(target_filenames[0])
        num_source_channels = self.num_input_channels
        num_target_channels = len(self.out_channel_names)

        for ch_iter in range(num_target_channels):

            for subj_iter in range(num_subjects):
                trg_img = nib.load(target_filenames[ch_iter][subj_iter])
                trg_img_data = trg_img.get_data()
                trg_img_data = self.preprocess_image(trg_img_data, self.out_channel_names[ch_iter])

                if self.use_slices == True:
                    if self.orientation is None:
                        print('orientation of slices is not specified')
                        return
                    elif self.orientation == 'coronal':
                        slice_idxs = np.asarray(range(0, trg_img_data.shape[2]))

                    trg_img_data = self.map_labels(trg_img_data, self.labels)
                    curr_channel = self.out_channel_names[ch_iter]
                    self.val_trg_image_dict[curr_channel].append(trg_img_data)
                    self.val_trg_slice_indices_dict[curr_channel].append(slice_idxs)
                    print(trg_img_data.shape)
                    print(slice_idxs.shape)

                else:

                    fg_idxs = np.argwhere(trg_img_data > 0)
                    trg_fg = trg_img_data[trg_img_data > 0].flatten()
                    trg_fg = np.reshape(trg_fg, (trg_fg.shape[0], 1))

                    print(trg_img_data.shape)
                    print(fg_idxs.shape)
                    curr_channel = self.out_channel_names[ch_iter]
                    self.val_trg_image_dict[curr_channel].append(trg_img_data)
                    self.val_trg_fg_indices_dict[curr_channel].append(fg_idxs)
                    self.val_trg_fg_labels_dict[curr_channel].append(trg_fg)

        for ch_iter in range(num_source_channels):

            curr_channel = self.channel_names[ch_iter]

            for subj_iter in range(num_subjects):
                src_img = nib.load(source_filenames[ch_iter][subj_iter])
                src_img_data = src_img.get_data()
                src_img_data = self.preprocess_image(src_img_data, channel_name=curr_channel)
                self.val_src_image_dict[curr_channel].append(src_img_data)

    def preprocess_image(self, in_img_data, channel_name):

        pad_width = ((self.feature_shape[0], self.feature_shape[0]),
                     (self.feature_shape[1], self.feature_shape[1]),
                     (self.feature_shape[2], self.feature_shape[2]))
        if self.use_slices == False:
            in_img_data = np.pad(in_img_data, pad_width, 'constant')

        if 'seg' in channel_name:
            # in_img_data = self.map_labels(in_img_data,self.labels )
            print('processing label image, do not do anything')

        else:

            in_img_data = in_img_data.astype(float)

            if 't1w' in channel_name:
                # it is a t1w image.
                if self.wmp_standardize == True:
                    if self.rob_standardize == True:
                        print('robnorm True wmp true')

                        in_img_data = robust_normalize(in_img_data)
                        in_img_data = wm_peak_normalize((in_img_data))

                    else:
                        print('robnorm False wmp true')
                        in_img_data = wm_peak_normalize((in_img_data))
                else:
                    print('robnorm true, wmp false')
                    in_img_data = robust_normalize(in_img_data)

                in_img_data[in_img_data > 255] = 255
                in_img_data = in_img_data / 255

            elif 't2w' in channel_name:
                # if self.wmp_standardize == True:
                #     in_img_data = intensity_standardize_utils.wm_peak_normalize_t2w(in_img_data)
                # else:
                in_img_data = wm_peak_normalize_t2w(in_img_data)

                # in_img_data[in_img_data > 255] = 255
                # in_img_data = in_img_data / 255
                # in_img_data[in_img_data > 1.2] = 1.2
            elif 'tal' in channel_name:
                # this channel has atlas coordinates
                print('channel name ' + channel_name)
                in_img_data = in_img_data / 127
            else:
                in_img_data = in_img_data
        return in_img_data

    def create_rare_training_feature_array(self, src_filenames, src_array_name, trg_filenames, trg_array_name,
                                           indices_list=None, step_size=[4, 4, 4], is_src_label_img=False,
                                           is_trg_label_img=True):
        nb_features_per_subject = 1000000
        nb_subjects = len(src_filenames)
        nb_src_modalities = len(src_filenames[0])
        print(src_filenames[0])
        tmpimg = nib.load(src_filenames[0])
        # tmpseg = nib.load(seg_filenames[0])

        if is_src_label_img == True:
            image_dtype = tmpimg.get_data_dtype()
        else:
            image_dtype = np.dtype(np.float32)

        src_feature_array = self.data_storage.create_earray(self.data_storage.root, src_array_name,
                                                            tables.Atom.from_dtype(image_dtype),
                                                            shape=(0,) + self.feature_shape + (1,),
                                                            expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
        trg_feature_array = self.data_storage.create_earray(self.data_storage.root, trg_array_name,
                                                            tables.Atom.from_dtype(image_dtype),
                                                            shape=(0,) + self.feature_shape + (1,),
                                                            expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        if self.use_patches == True:
            src_index_array = self.data_storage.create_earray(self.data_storage.root, src_array_name + '_index',
                                                              tables.Int16Atom(), shape=(0, 3),
                                                              expectedrows=np.prod(
                                                                  nb_features_per_subject) * nb_subjects)

            trg_index_array = self.data_storage.create_earray(self.data_storage.root, trg_array_name + '_index',
                                                              tables.Int16Atom(), shape=(0, 3),
                                                              expectedrows=np.prod(
                                                                  nb_features_per_subject) * nb_subjects)

            for src_file, trg_file in zip(src_filenames, trg_filenames):
                (trg_features, trg_indices) = self.extract_training_patches(trg_file, intensity_threshold=0,
                                                                            step_size=step_size, indices=None,
                                                                            is_label_img=is_trg_label_img)

                rare_idxs = np.where((trg_features == 4) | (trg_features == 11) | (trg_features == 12) |
                                     (trg_features == 19) | (trg_features == 24) | (trg_features == 38))
                rare_patch_idx = np.unique(rare_idxs[0])
                rare_trg_features = trg_features[rare_patch_idx]
                rare_trg_indices = trg_indices[rare_patch_idx]

                (src_features, indices) = self.extract_training_patches(src_file, intensity_threshold=0,
                                                                        step_size=step_size, indices=rare_trg_indices,
                                                                        is_label_img=is_src_label_img)

                print(src_file + " src rare feature size ")
                print(src_features.shape)
                print(trg_file + " trg rare feature size ")
                print(rare_trg_features.shape)

                src_feature_array.append(src_features)
                src_index_array.append(indices)
                trg_feature_array.append(rare_trg_features)
                trg_index_array.append(indices)

                if indices_list == None:
                    indices_list = list()

                indices_list.append(indices)

        return src_feature_array, src_index_array, trg_feature_array, trg_index_array, indices_list

    def create_training_feature_array(self, image_filenames, array_name, indices_list, step_size, is_label_img,
                                      orientation):
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
                                                        shape=(0,) + self.feature_shape[0:-1] + (1,),
                                                        expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        if self.use_patches == True:
            index_array = self.data_storage.create_earray(self.data_storage.root, array_name + '_index',
                                                          tables.Int16Atom(), shape=(0, 3),
                                                          expectedrows=np.prod(nb_features_per_subject) * nb_subjects)
            if indices_list == None:
                print("No indices_list found")
                indices_list = list()
                for input_file in image_filenames:
                    (features, indices) = self.extract_training_patches(input_file,
                                                                        intensity_threshold=0,
                                                                        step_size=step_size,
                                                                        indices=None,
                                                                        is_label_img=is_label_img,
                                                                        channel_name=array_name)
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
                    (features, indices) = self.extract_training_patches(input_file,
                                                                        intensity_threshold=0,
                                                                        step_size=step_size,
                                                                        indices=curr_indices,
                                                                        is_label_img=is_label_img,
                                                                        channel_name=array_name)
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
            print(indices_list)
            if indices_list == None:
                indices_list = list()
                if self.dim == 2:
                    for input_file in image_filenames:
                        (features, indices) = self.extract_training_slices(input_file, intensity_threshold=0,
                                                                           indices=None, is_label_img=is_label_img,
                                                                           orientation=orientation)
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
                                                                       indices=curr_indices, is_label_img=is_label_img,
                                                                       orientation=orientation)
                    print("indices shape is ")
                    print(indices.shape)
                    feature_array.append(features)
                    index_array.append(curr_indices)
                    print(input_file + " features extract size ")
                    print(features.__class__)

        return feature_array, index_array, indices_list

    def create_training_label_array(self, target_label_list, array_name, indices_list):
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
            curr_label_array = np.ones(curr_indices.shape) * curr_label
            label_array.append(curr_label_array)
            index_array.append(curr_indices)

    def extract_training_slices(self, in_img_file, intensity_threshold, indices, is_label_img, orientation):
        if indices is not None:
            (slices, indices) = self.extract_slices(in_img_file, intensity_threshold, is_label_img=is_label_img,
                                                    indices=indices, orientation=orientation)
            return slices, indices
        else:
            (slices, indices) = self.extract_slices(in_img_file, intensity_threshold, is_label_img=is_label_img,
                                                    indices=indices, orientation=orientation)
            training_slices = slices
            training_indices = indices
            return training_slices, np.int32(training_indices)

    def extract_slices(self, in_img_file, intensity_threshold, indices=None, is_label_img=False, orientation='coronal'):
        in_img = nib.load(in_img_file)

        # resize to 256x
        # white matter peak set to 200 and divide by 255

        if is_label_img == False:
            in_img_data = in_img.get_data().astype(float)
            if self.wmp_standardize == True:
                in_img_data = wm_peak_normalize(in_img_data)
                in_img_data[in_img_data > 255] = 255
                in_img_data = in_img_data / 255
            else:
                in_img_data = robust_normalize(in_img_data)
                in_img_data = in_img_data / 255
            print("Image max is :" + str(in_img_data.max()))
        else:
            in_img_data = in_img.get_data()

        if orientation == 'axial':
            in_img_data = np.transpose(in_img_data, (0, 2, 1))
        elif orientation == 'sagittal':
            in_img_data = np.transpose(in_img_data, (2, 1, 0))

        if indices is not None:
            slices = []
            for z_index in indices:
                slices.append(in_img_data[:, :, z_index])
        else:
            slices = []
            indices = range(in_img_data.shape[2])
            print(in_img_data.shape[2])
            for z_index in indices:
                print(z_index)
                slices.append(in_img_data[:, :, z_index])
            slices = np.asarray(slices)
            indices = np.asarray(indices)
            # add channel as a dimension for keras
            newshape = list(slices.shape)
            newshape.append(1)
            slices = np.reshape(slices, newshape)
            indices = indices.reshape(-1, 1)

        return slices, indices

    def extract_training_patches(self, in_img_file, intensity_threshold, step_size, indices, is_label_img,
                                 channel_name):
        if indices is not None:
            print('channel name is ' + channel_name)
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                         is_label_img=is_label_img, indices=indices,
                                                         channel_name=channel_name)
            return patches, indices
        else:
            (patches, indices, _) = self.extract_patches(in_img_file, intensity_threshold, step_size,
                                                         is_label_img=is_label_img, indices=indices,
                                                         channel_name=channel_name)
            # (seg_patches, seg_indices, _) = self.extract_patches(seg_img_file, intensity_threshold, step_size,
            #                                                   is_label_img=True, indices=indices)
            training_patches = patches
            training_indices = indices
            return training_patches, np.int32(training_indices)

    def extract_patches(self, in_img_file, intensity_threshold, step_size, is_label_img=False,
                        indices=None, robnorm=None, channel_name=None):
        # pad the images by patch shape
        in_img = nib.load(in_img_file)
        print('in extract_patches channel_name is ' + channel_name)

        # white matter peak set to 200 and divide by 255
        if is_label_img == False:
            in_img_data = in_img.get_data().astype(float)

            if 't1w' in channel_name:
                print('channel name' + channel_name)
                # it is a t1w image.
                if self.wmp_standardize == True:
                    if self.rob_standardize == False:
                        print('robnorm False')
                        # in_img_data = robust_normalize(in_img_data)
                        in_img_data = wm_peak_normalize(in_img_data)
                    else:
                        print('robnorm True')
                        in_img_data = robust_normalize(in_img_data)
                        in_img_data = wm_peak_normalize(in_img_data)
                else:
                    if self.rob_standardize == False:
                        print('robnorm False wmp False')
                    else:
                        print('robnorm True wmp False')
                        in_img_data = robust_normalize(in_img_data)
                    in_img_data[in_img_data > 255] = 255
                    in_img_data = in_img_data / 255

            elif 'tal' in channel_name:
                # this channel has atlas coordinates
                print('channel name ' + channel_name)
                in_img_data = in_img_data / 127
            else:
                in_img_data = in_img_data

        else:
            in_img_data = in_img.get_data()
            in_img_data = self.map_labels(in_img_data, self.labels)

        padding0 = (self.feature_shape[0] + step_size[0] + 1, self.feature_shape[0] + step_size[0] + 1)
        padding1 = (self.feature_shape[1] + step_size[1] + 1, self.feature_shape[1] + step_size[1] + 1)
        padding2 = (self.feature_shape[2] + step_size[2] + 1, self.feature_shape[2] + step_size[2] + 1)

        in_img_data_pad = np.pad(in_img_data, (padding0, padding1, padding2), 'constant', constant_values=0)
        padded_img_size = in_img_data_pad.shape

        if indices is not None:

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
        if (is_label_img == False) and (self.preprocessing == True):
            print('Unit normalizing')
            orig_shape = patches.shape
            patches = patches.reshape(orig_shape[0], -1)
            patches = sklearn.preprocessing.scale(patches)
            patches = patches.reshape(orig_shape)

        # add channel as a dimension for keras
        newshape = list(patches.shape)
        newshape.append(1)
        patches = np.reshape(patches, newshape)

        if indices is None:
            chdim = len(patches.shape) - 1
            psum = np.sum(patches, axis=chdim)
            for a in range(chdim - 1):
                psum = np.sum(psum, axis=chdim - a - 1)

            nz_idxs = np.where(psum != 0)[0]
            print(str(patches.shape[0]) + ' total patches removed')
            patches = patches[nz_idxs, :]
            print(str(patches.shape[0]) + ' after patches removed')

            indices = np.concatenate((idx_x.reshape(-1, 1), idx_y.reshape(-1, 1), idx_z.reshape(-1, 1)), axis=1)
            indices = indices[nz_idxs, :]

        print('Patches shape is ' + str(patches.shape))
        print('Indices shape is ' + str(indices.shape))
        return patches, np.int32(indices), padded_img_size

    def build_image_from_patches(self, in_patches, indices, padded_img_size, patch_crop_size, step_size):
        """
        Note: patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]

        """
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
                np.multiply(np.reshape(in_patches[patch_iter, :], out_feature_shape), patch_mask)

            count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2]] += patch_mask

        out_img_data = np.divide(out_img_data, count_img_data)
        out_img_data[np.isnan(out_img_data)] = 0
        # remove the padding
        unpadded_img_size = padded_img_size - np.multiply(np.asarray(self.feature_shape[:-1]) + step_size + 1, 2)
        padding = np.asarray(self.feature_shape[:-1]) + step_size + 1

        print("padding is " + str(np.multiply(np.asarray(self.feature_shape[:-1]) + step_size + 1, 2)))
        print("unpadded image size is " + str(unpadded_img_size))

        out_img_data = out_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                       padding[1]:padding[1] + unpadded_img_size[1],
                       padding[2]:padding[2] + unpadded_img_size[2]]
        count_img_data = count_img_data[padding[0]:padding[0] + unpadded_img_size[0],
                         padding[1]:padding[1] + unpadded_img_size[1],
                         padding[2]:padding[2] + unpadded_img_size[2]]
        return out_img_data, count_img_data

    def build_seg_from_patches(self, in_patches, indices, padded_img_size, patch_crop_size, step_size):
        """
        Note: patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]

        """
        out_img_data = np.zeros(padded_img_size)
        count_img_data = np.zeros(padded_img_size)

        out_feature_shape = in_patches.shape[1:5]

        idx_x = indices[:, 0]
        idx_y = indices[:, 1]
        idx_z = indices[:, 2]

        patch_mask = np.zeros(out_feature_shape)
        patch_mask[0 + patch_crop_size[0]: out_feature_shape[0] - patch_crop_size[0],
        0 + patch_crop_size[1]: out_feature_shape[1] - patch_crop_size[1],
        0 + patch_crop_size[2]: out_feature_shape[2] - patch_crop_size[2], :] = 1

        for patch_iter in range(len(idx_x)):
            out_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2], :] += \
                np.multiply(np.reshape(in_patches[patch_iter, :], out_feature_shape), patch_mask)

            count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + out_feature_shape[0],
            idx_y[patch_iter]:idx_y[patch_iter] + out_feature_shape[1],
            idx_z[patch_iter]:idx_z[patch_iter] + out_feature_shape[2], :] += patch_mask

        out_img_data = np.divide(out_img_data, count_img_data)
        out_img_data[np.isnan(out_img_data)] = 0

        # remove the padding
        # unpadded_img_size = padded_img_size[0:3] - np.multiply(self.feature_shape, 2)

        unpadded_img_size = padded_img_size[0:3] - np.multiply(np.asarray(self.feature_shape[:-1]) + step_size + 1, 2)
        padding = np.asarray(self.feature_shape[:-1]) + step_size + 1

        print("padding is " + str(np.multiply(np.asarray(self.feature_shape[:-1]) + step_size + 1, 2)))
        print("unpadded image size is " + str(unpadded_img_size))

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

    def build_seg_from_slices(self, in_slices, indices):
        out_img_data = np.zeros(in_slices.shape[1:])

    def generate_src_trg_training_data(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                                       target_label_list=None, step_size=None, ):
        print("number of channels is " + str(self.num_input_channels))
        for iter_channel in range(self.num_input_channels):
            source_filenames_channel = source_filenames[iter_channel]
            print(source_filenames_channel)
            print('Creating source image patches.')

            if self.use_patches == True:
                if iter_channel == 0:
                    (_, _, indices_list) = self.create_training_feature_array(source_filenames_channel,
                                                                              'src' + self.channel_names[iter_channel],
                                                                              indices_list=None, step_size=step_size,
                                                                              is_label_img=is_src_label_img)
                else:
                    (_, _, _) = self.create_training_feature_array(source_filenames_channel,
                                                                   'src' + str(self.channel_names[iter_channel]),
                                                                   indices_list=indices_list, step_size=step_size,
                                                                   is_label_img=is_src_label_img, )
            else:
                (_, _, indices_list) = self.create_training_feature_array(source_filenames, 'src',
                                                                          indices_list=None, step_size=None,
                                                                          is_label_img=is_src_label_img)
                self.create_training_label_array(target_label_list, 'trg', indices_list)

        print('NCreating target image patches.')
        print("number of output channels is " + str(len(self.out_channel_names)))
        for iter_channel in range(len(self.out_channel_names)):
            target_filenames_channel = target_filenames[iter_channel]
            self.create_training_feature_array(target_filenames_channel, 'trg' + self.out_channel_names[iter_channel],
                                               indices_list, step_size=step_size, is_label_img=is_trg_label_img)

            # print('Creating rare source and target  patches')
            # (_, _, _, _, rare_indices_list) = self.create_rare_training_feature_array(source_filenames,
            #                                                                      'src_rare',
            #                                                                      target_filenames,
            #                                                                      'trg_rare',
            #                                                                      indices_list=None,
            #                                                                      step_size=step_size,
            #                                                                      is_src_label_img=is_src_label_img,
            #                                                                      is_trg_label_img=is_trg_label_img)

    def generate_src_trg_validation_data(self, source_filenames, target_filenames, is_src_label_img, is_trg_label_img,
                                         target_label_list=None, step_size=None):

        print("number of channels is " + str(self.num_input_channels))
        for iter_channel in range(self.num_input_channels):
            source_filenames_channel = source_filenames[iter_channel]
            print(source_filenames_channel)

            print('Creating source image patches.')
            if self.use_patches == True:
                if iter_channel == 0:
                    (_, _, indices_list) = self.create_training_feature_array(source_filenames_channel,
                                                                              'val_src' + str(
                                                                                  self.channel_names[iter_channel]),
                                                                              indices_list=None, step_size=step_size,
                                                                              is_label_img=is_src_label_img, )
                else:
                    (_, _, _) = self.create_training_feature_array(source_filenames_channel,
                                                                   'val_src' + str(self.channel_names[iter_channel]),
                                                                   indices_list=indices_list, step_size=step_size,
                                                                   is_label_img=is_src_label_img, )
            else:
                (_, _, indices_list) = self.create_training_feature_array(source_filenames, 'valsrc',
                                                                          indices_list=None, step_size=None,
                                                                          is_label_img=is_src_label_img)
                self.create_training_label_array(target_label_list, 'valtrg', indices_list)

        print("number of output channels is " + str(len(self.out_channel_names)))
        for iter_channel in range(len(self.out_channel_names)):
            target_filenames_channel = target_filenames[iter_channel]
            self.create_training_feature_array(target_filenames_channel,
                                               'val_trg' + self.out_channel_names[iter_channel],
                                               indices_list, step_size=step_size, is_label_img=is_trg_label_img)

    def generate_src_trg_training_slice_data(self, source_filenames, target_filenames, is_src_label_img,
                                             is_trg_label_img, orientation):

        print("number of channels is " + str(self.num_input_channels))
        for iter_channel in range(self.num_input_channels):
            source_filenames_channel = source_filenames[iter_channel]
            print(source_filenames_channel)

            print('Creating source image patches.')
            if iter_channel == 0:
                (_, _, indices_list) = self.create_training_feature_array(source_filenames_channel,
                                                                          'src' + self.channel_names[iter_channel],
                                                                          indices_list=None, step_size=None,
                                                                          is_label_img=is_src_label_img,
                                                                          orientation=orientation)
            else:
                (_, _, _) = self.create_training_feature_array(source_filenames_channel,
                                                               'src' + str(self.channel_names[iter_channel]),
                                                               indices_list=indices_list, step_size=None,
                                                               is_label_img=is_src_label_img, orientation=orientation)
        print('NCreating target image patches.')
        print("number of output channels is " + str(len(self.out_channel_names)))
        for iter_channel in range(len(self.out_channel_names)):
            target_filenames_channel = target_filenames[iter_channel]
            self.create_training_feature_array(target_filenames_channel, 'trg' + self.out_channel_names[iter_channel],
                                               indices_list, step_size=None, is_label_img=is_trg_label_img,
                                               orientation=orientation)

    def generate_src_trg_validation_slice_data(self, source_filenames, target_filenames, is_src_label_img,
                                               is_trg_label_img, orientation):

        print("number of channels is " + str(self.num_input_channels))
        for iter_channel in range(self.num_input_channels):
            source_filenames_channel = source_filenames[iter_channel]
            print(source_filenames_channel)

            print('Creating source image patches.')
            if iter_channel == 0:
                (_, _, indices_list) = self.create_training_feature_array(source_filenames_channel,
                                                                          'val_src' + self.channel_names[iter_channel],
                                                                          indices_list=None, step_size=None,
                                                                          is_label_img=is_src_label_img,
                                                                          orientation=orientation
                                                                          )
            else:
                (_, _, _) = self.create_training_feature_array(source_filenames_channel,
                                                               'val_src' + str(self.channel_names[iter_channel]),
                                                               indices_list=indices_list, step_size=None,
                                                               is_label_img=is_src_label_img,
                                                               orientation=orientation)

        print('NCreating target image patches.')
        print("number of output channels is " + str(len(self.out_channel_names)))
        for iter_channel in range(len(self.out_channel_names)):
            target_filenames_channel = target_filenames[iter_channel]
            self.create_training_feature_array(target_filenames_channel,
                                               'val_trg' + self.out_channel_names[iter_channel],
                                               indices_list, step_size=None, is_label_img=is_trg_label_img,
                                               orientation=orientation)

    def extract_training_image_slices(self, source_filename, target_filename, intensity_threshold=0, step_size=None,
                                      dataset=None):

        imgs = nib.load(source_filename)
        vecs = nib.load(target_filename)

        img_data = imgs.get_data()
        img_data = img_data[:, :, 0]

        vec_data = vecs.get_data()
        vec_data = vec_data[:, :, 0, :]

        padding0 = (self.feature_shape[0] + step_size[0] + 1, self.feature_shape[0] + step_size[0] + 1)
        padding1 = (self.feature_shape[1] + step_size[1] + 1, self.feature_shape[1] + step_size[1] + 1)
        padding2 = (0, 0)
        img_data_pad = np.pad(img_data, (padding0, padding1), 'constant', constant_values=0)
        vec_data_pad = np.pad(vec_data, (padding0, padding1, padding2), 'constant', constant_values=0)

        padded_img_size = img_data_pad.shape
        intensity_threshold = 0

        (idx_x_fg, idx_y_fg, idx_z_fg) = np.where(vec_data_pad > intensity_threshold)
        min_idx_x_fg = np.min(idx_x_fg) - step_size[0]
        max_idx_x_fg = np.max(idx_x_fg) + step_size[0]
        min_idx_y_fg = np.min(idx_y_fg) - step_size[1]
        max_idx_y_fg = np.max(idx_y_fg) + step_size[1]

        sampled_x = np.arange(min_idx_x_fg, max_idx_x_fg, step_size[0])
        sampled_y = np.arange(min_idx_y_fg, max_idx_y_fg, step_size[1])

        idx_x, idx_y = np.meshgrid(sampled_x, sampled_y, sparse=False, indexing='ij')
        idx_x = idx_x.flatten()
        idx_y = idx_y.flatten()
        # idx_z = idx_z.flatten()

        patches = []
        trg_patches = []
        indices_x = []
        indices_y = []
        # take half the patches (because we have no memory)

        print('dataset is for ' + dataset)
        if dataset == 'train':
            iterrange = range(0, len(idx_x) // 3)
        elif dataset == 'val':
            iterrange = range((len(idx_x) // 3 + 1), (len(idx_x) // 3) + (len(idx_x) // 5))
        print('iter range is ' + str(iterrange[0]) + '_' + str(iterrange[-1]))

        for patch_iter in iterrange:
            curr_patch = img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                         idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1]]
            vec_patch = vec_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + self.feature_shape[0],
                        idx_y[patch_iter]:idx_y[patch_iter] + self.feature_shape[1], 0:2]
            if vec_patch.mean() != 0:
                print(patch_iter * 100.0 / len(idx_x))
                patches.append(curr_patch)
                trg_patches.append(vec_patch)
                indices_x.append(idx_x[patch_iter])
                indices_y.append(idx_y[patch_iter])
                # save the patches in out.mgz files
                # currshape = list(curr_patch.shape)
                # currshape.append(1)
                # vecshape = list(vec_patch.shape[0:2])
                # vecshape.append(1)
                # vecshape.append(2)

        patches = np.asarray(patches)
        newshape = list(patches.shape)
        newshape.append(1)
        indices = np.vstack((np.asarray(indices_x), np.asarray(indices_y))).T
        patches = np.reshape(patches, newshape)
        print(indices.shape)
        return patches, np.asarray(trg_patches), indices

    def generate_src_trg_slice_training_validation(self, source_filename, target_filename, step_size=None):
        nb_features_per_subject = 10000
        nb_subjects = 50000
        tmpimg = nib.load(source_filename)
        # tmpseg = nib.load(seg_filenames[0])
        image_dtype = np.dtype(np.float32)

        feature_array = self.data_storage.create_earray(self.data_storage.root, 'src',
                                                        tables.Atom.from_dtype(image_dtype),
                                                        shape=(0,) + self.feature_shape + (1,),
                                                        expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        target_array = self.data_storage.create_earray(self.data_storage.root, 'trg',
                                                       tables.Atom.from_dtype(image_dtype),
                                                       shape=(0,) + self.feature_shape + (self.num_outputs,),
                                                       expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        index_array = self.data_storage.create_earray(self.data_storage.root, 'src' + '_index',
                                                      tables.Int16Atom(), shape=(0, 2),
                                                      expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        val_feature_array = self.data_storage.create_earray(self.data_storage.root, 'src_validation',
                                                            tables.Atom.from_dtype(image_dtype),
                                                            shape=(0,) + self.feature_shape + (1,),
                                                            expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        val_target_array = self.data_storage.create_earray(self.data_storage.root, 'trg_validation',
                                                           tables.Atom.from_dtype(image_dtype),
                                                           shape=(0,) + self.feature_shape + (self.num_outputs,),
                                                           expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        val_index_array = self.data_storage.create_earray(self.data_storage.root, 'val' + '_index',
                                                          tables.Int16Atom(), shape=(0, 2),
                                                          expectedrows=np.prod(nb_features_per_subject) * nb_subjects)

        (slices, slice_labels, indices) = self.extract_training_image_slices(source_filename,
                                                                             target_filename,
                                                                             intensity_threshold=0,
                                                                             step_size=step_size, dataset='train')
        feature_array.append(slices)
        target_array.append(slice_labels)
        index_array.append(indices)

        (val_slices, val_slice_labels, val_indices) = self.extract_training_image_slices(source_filename,
                                                                                         target_filename,
                                                                                         intensity_threshold=0,
                                                                                         step_size=step_size,
                                                                                         dataset='val')
        print('val slice shape ' + str(val_slices.shape))
        print('val target shape ' + str(val_slice_labels.shape))
        print('val idx shape ' + str(val_indices.shape))

        val_feature_array.append(val_slices)
        val_target_array.append(val_slice_labels)
        val_index_array.append(val_indices)

        return feature_array, index_array, indices

    def dynamic_seg_training_generator_singlechannel(self, batch_size):
        while True:
            x_list = list()
            y_list = list()
            trg_channel_name = self.out_channel_names[0]
            img_index_list = range(len(self.trg_image_dict[trg_channel_name]))
            #   list(range(self.data_storage.root.trgseg.shape[0]))

            # Randomly generate batch_size of subject_ids.
            num_subjects = len(img_index_list)
            subj_idxs = np.random.choice(range(0, num_subjects), batch_size)
            batch_index_list = list()
            for subj_idx in subj_idxs:
                p_idxs = random.sample(xrange(len(self.trg_fg_indices_dict[trg_channel_name][subj_idx])), 1)
                for p_idx in p_idxs:
                    [px, py, pz] = self.trg_fg_indices_dict[trg_channel_name][subj_idx][p_idx]
                    batch_index_list.append([subj_idx, px, py, pz])

            batch_index_array = np.array(batch_index_list)
            num_orig_mprage = (batch_size)

            orig_mprage_idxs = range(0, num_orig_mprage)

            orig_mprage_batch_idxs = batch_index_array[orig_mprage_idxs]

            for b_iter in range(orig_mprage_batch_idxs.shape[0]):
                subj_iter = orig_mprage_batch_idxs[b_iter][0]
                px = orig_mprage_batch_idxs[b_iter][1]
                py = orig_mprage_batch_idxs[b_iter][2]
                pz = orig_mprage_batch_idxs[b_iter][3]
                x_t1w_patch = self.src_image_dict['t1w'][subj_iter][
                              int(px) - self.feature_shape[0] / 2: int(px) + self.feature_shape[0] / 2,
                              int(py) - self.feature_shape[1] / 2: int(py) + self.feature_shape[1] / 2,
                              int(pz) - self.feature_shape[2] / 2:int(pz) + self.feature_shape[2] / 2]

                y_patch = self.trg_image_dict['seg'][subj_iter][
                          int(px) - self.feature_shape[0] / 2: int(px) + self.feature_shape[0] / 2,
                          int(py) - self.feature_shape[1] / 2: int(py) + self.feature_shape[1] / 2,
                          int(pz) - self.feature_shape[2] / 2:int(pz) + self.feature_shape[2] / 2]
                x_list.append(x_t1w_patch)
                y_list.append(y_patch)

            x_array = np.stack(x_list, axis=0)
            y_array = np.stack(y_list, axis=0)
            x_array = np.reshape(x_array, x_array.shape + (1,))
            y_array = np.reshape(y_array, y_array.shape + (1,))
            x, y = self.convert_data(x_array, y_array, self.n_labels, self.labels)

            yield x, y

    def dynamic_seg_validation_generator_singlechannel(self, batch_size):
        while True:
            x_list = list()
            y_list = list()
            trg_channel_name = self.out_channel_names[0]
            img_index_list = range(len(self.val_trg_image_dict[trg_channel_name]))
            #   list(range(self.data_storage.root.trgseg.shape[0]))

            # Randomly generate batch_size of subject_ids.
            num_subjects = len(img_index_list)
            subj_idxs = np.random.choice(range(0, num_subjects), batch_size)
            batch_index_list = list()
            for subj_idx in subj_idxs:
                p_idxs = random.sample(xrange(len(self.val_trg_fg_indices_dict[trg_channel_name][subj_idx])), 1)
                for p_idx in p_idxs:
                    [px, py, pz] = self.val_trg_fg_indices_dict[trg_channel_name][subj_idx][p_idx]
                    batch_index_list.append([subj_idx, px, py, pz])

            for b_index in batch_index_list:
                subj_iter = b_index[0]
                px = b_index[1]
                py = b_index[2]
                pz = b_index[3]
                x_t1w_patch = self.val_src_image_dict['t1w'][subj_iter][
                              int(px) - self.feature_shape[0] / 2: int(px) + self.feature_shape[0] / 2,
                              int(py) - self.feature_shape[1] / 2: int(py) + self.feature_shape[1] / 2,
                              int(pz) - self.feature_shape[2] / 2:int(pz) + self.feature_shape[2] / 2]

                y_patch = self.val_trg_image_dict['seg'][subj_iter][
                          int(px) - self.feature_shape[0] / 2: int(px) + self.feature_shape[0] / 2,
                          int(py) - self.feature_shape[1] / 2: int(py) + self.feature_shape[1] / 2,
                          int(pz) - self.feature_shape[2] / 2:int(pz) + self.feature_shape[2] / 2]
                x_list.append(x_t1w_patch)
                y_list.append(y_patch)

            x_array = np.stack(x_list, axis=0)
            y_array = np.stack(y_list, axis=0)
            x_array = np.reshape(x_array, x_array.shape + (1,))
            y_array = np.reshape(y_array, y_array.shape + (1,))
            x, y = self.convert_data(x_array, y_array, self.n_labels, self.labels)

            yield x, y

    def training_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.srct1w[index])
                y_list.append(self.data_storage.root.trg[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list

    def validation_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.val_srct1w[index])
                y_list.append(self.data_storage.root.val_trg[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list

    def training_generator_t1beta(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.srct1w[index])
                y_list.append(self.data_storage.root.trgt1beta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            # normalize t1
            y_list = y_list / 5000
            yield x_list, y_list

    def validation_generator_t1beta(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.val_srct1w[index])
                y_list.append(self.data_storage.root.val_trgt1beta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)

            y_list = y_list / 5000
            yield x_list, y_list

    def training_generator_t2beta(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.srct1w[index])
                y_list.append(self.data_storage.root.trgt2beta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            # normalize t1
            y_list = y_list / 3000
            yield x_list, y_list

    def validation_generator_t2beta(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.val_srct1w[index])
                y_list.append(self.data_storage.root.val_trgt2beta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            y_list = y_list / 3000
            yield x_list, y_list

    def training_generator_pdbeta(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.srct1w[index])
                y_list.append(self.data_storage.root.trgpdbeta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            y_list = y_list / 15000
            yield x_list, y_list

    def validation_generator_pdbeta(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.val_srct1w[index])
                y_list.append(self.data_storage.root.val_trgpdbeta[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)

            y_list = y_list / 15000
            yield x_list, y_list

    def map_labels(self, input_label_patch, input_label_list):
        output_label_patch = np.zeros(input_label_patch.shape)
        # 0th label is always 0
        for out_label, in_label in enumerate(input_label_list):
            output_label_patch[input_label_patch == in_label] = out_label

        return output_label_patch

    def map_inv_labels(self, input_label_patch, input_label_list):
        output_label_patch = np.zeros(input_label_patch.shape)
        # 0th label is always 0
        for curr_label in range(len(input_label_list)):
            output_label_patch[input_label_patch == curr_label] = input_label_list[curr_label]

        return output_label_patch

    def get_multi_class_labels(self, data, n_labels):
        """Translates a label map into a set of binary labels.

        Args:
            data: numpy array containing the label map with shape: (n_samples, 1, ...).
            n_labels: number of labels.
            labels: integer values of the labels.

        Returns:
            A binary numpy array of shape: (n_samples, n_labels, ...)

        """
        labels = range(n_labels)

        new_shape = list(data.shape[0:-1]) + [n_labels]  # [data.shape[0], n_labels] + list(data.shape[2:])
        y = np.zeros(new_shape, np.int8)
        if len(data.shape) == 4:  # 2D data
            for label_index in range(n_labels):
                if labels is not None:
                    y[:, :, :, label_index][data[:, :, :, 0] == labels[label_index]] = 1
                else:
                    y[:, :, :, label_index][data[:, :, :, 0] == (label_index + 1)] = 1
            return y
        elif len(data.shape) == 5:  # 3D data
            for label_index in range(n_labels):
                if labels is not None:
                    y[:, :, :, :, label_index][data[:, :, :, :, 0] == labels[label_index]] = 1
                else:
                    y[:, :, :, :, label_index][data[:, :, :, :, 0] == (label_index + 1)] = 1
            return y

    def convert_data(self, x, y, n_labels=1, labels=None):
        if n_labels == 1:
            y[y > 0] = 1
        elif n_labels > 1:
            y = self.get_multi_class_labels(y, n_labels=n_labels)
        return x, y

    def seg_training_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))

        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            random.shuffle(rare_index_list)
            for index in index_list[:num_remaining_patches]:
                x_list.append(self.data_storage.root.src[index])
                y_list.append(self.data_storage.root.trg[index])

            for index in rare_index_list[:num_rare_patches]:
                x_list.append(self.data_storage.root.src_rare[index])
                y_list.append(self.data_storage.root.trg_rare[index])

            x = np.asarray(x_list)
            y = np.asarray(y_list)
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)

            yield x, y

    def seg_training_generator_multichannel(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            # random.shuffle(rare_index_list)
            xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:batch_size], :],
                              self.data_storage.root.srctalx[index_list[0:batch_size], :],
                              self.data_storage.root.srctaly[index_list[0:batch_size], :],
                              self.data_storage.root.srctalz[index_list[0:batch_size], :]), axis=-2)
            xfeat = xfeat.reshape(xfeat.shape[:-1])
            y = self.data_storage.root.trg[index_list[0:batch_size], :]
            x = xfeat
            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            # y = self.map_inv_labels(y, self.labels)
            yield x, y

    def seg_validation_generator_multichannel(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            # random.shuffle(rare_index_list)
            xfeat = np.stack((self.data_storage.root.val_srct1w[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctalx[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctaly[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctalz[index_list[0:batch_size], :]), axis=-2)
            xfeat = xfeat.reshape(xfeat.shape[:-1])
            y = self.data_storage.root.val_trg[index_list[0:batch_size], :]
            x = xfeat
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_training_generator_singlechannel(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            # random.shuffle(rare_index_list)
            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctalx[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctaly[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctalz[index_list[0:batch_size],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])
            x = self.data_storage.root.srct1w[index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]
            # x = xfeat
            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            # y = self.map_inv_labels(y, self.labels)
            yield x, y

    def seg_validation_generator_singlechannel(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            # random.shuffle(rare_index_list)
            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctalx[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctaly[index_list[0:batch_size],:],
            #                   self.data_storage.root.srctalz[index_list[0:batch_size],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])
            x = self.data_storage.root.val_srct1w[index_list[0:batch_size], :]
            y = self.data_storage.root.val_trgseg[index_list[0:batch_size], :]
            # x = xfeat
            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            # y = self.map_inv_labels(y, self.labels)
            yield x, y

    def apply_flash(self, theta_flash, batch_pd, batch_t1, batch_t2):
        batch_flash = np.zeros(batch_pd.shape)
        batch_flash[batch_pd > 0] = np.exp(np.log(batch_pd[batch_pd > 0]) + theta_flash[0] +
                                           (theta_flash[1] / batch_t1[batch_pd > 0]) +
                                           (theta_flash[2] / batch_t2[batch_pd > 0]))
        batch_flash = np.nan_to_num(batch_flash)
        batch_flash[batch_flash < 0] = 0.01
        batch_flash[batch_flash > 1] = 1
        batch_flash = scipy.signal.medfilt(batch_flash, kernel_size=[1, 3, 3, 3, 1])
        return batch_flash

    def apply_mprage(self, theta_mprage, batch_pd, batch_t1, batch_t2):
        batch_mprage = np.zeros(batch_pd.shape)
        batch_mprage[batch_pd > 0] = np.exp(np.log(batch_pd[batch_pd > 0]) + theta_mprage[0] +
                                            (theta_mprage[1] * batch_t1[batch_pd > 0]) +
                                            (theta_mprage[2] * batch_t1[batch_pd > 0] ** 2))
        batch_mprage = np.nan_to_num(batch_mprage)
        batch_mprage[batch_mprage < 0] = 0.01
        batch_mprage[batch_mprage > 1] = 1.0
        return batch_mprage

    def seg_training_generator_multichannel_nmr(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_orig_patches = batch_size // 4
        num_flash_patches = batch_size // 2
        num_mprage_patches = batch_size - (num_orig_patches + num_flash_patches)

        while True:
            x_list = list()
            y_list = list()
            # t= time.time()
            random.shuffle(index_list)
            # print(time.time() - t)

            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalx[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctaly[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalz[index_list[0:num_orig_patches],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t2 = time.time()
            x_t1w = self.data_storage.root.srct1w[index_list[0:batch_size], :]
            x_talx = self.data_storage.root.srctalx[index_list[0:batch_size], :]
            x_taly = self.data_storage.root.srctaly[index_list[0:batch_size], :]
            x_talz = self.data_storage.root.srctalz[index_list[0:batch_size], :]

            x_pd = self.data_storage.root.srcnmrpd[index_list[0:batch_size], :]
            x_t1 = self.data_storage.root.srcnmrt1[
                   index_list[0:batch_size], :]
            x_t2 = self.data_storage.root.srcnmrt2[
                   index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]

            # print(time.time() - t2)
            xfeat = np.stack((x_t1w[0:num_orig_patches, :],
                              x_talx[0:num_orig_patches, :],
                              x_taly[0:num_orig_patches, :],
                              x_talz[0:num_orig_patches, :]), axis=-2)
            xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t3 = time.time()
            xflash_pd = x_pd[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t1 = x_t1[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t2 = x_t2[num_orig_patches:num_orig_patches + num_flash_patches, :]

            flash_idx = np.random.randint(0, self.theta_flash.shape[0])
            curr_theta_flash = self.theta_flash[flash_idx, :]

            xflash = self.apply_flash(curr_theta_flash, xflash_pd, xflash_t1, xflash_t2)
            xflash_feat = np.stack((xflash,
                                    x_talx[num_orig_patches:num_orig_patches + num_flash_patches],
                                    x_taly[num_orig_patches:num_orig_patches + num_flash_patches],
                                    x_talz[num_orig_patches:num_orig_patches + num_flash_patches],
                                    ), axis=-2)
            xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # print(time.time() - t3)
            # xflash_feat = np.stack((xflash,
            #                   self.data_storage.root.srctalx[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctaly[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctalz[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :]), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # t4 = time.time()
            xmprage_pd = x_pd[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t1 = x_t1[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t2 = x_t2[num_orig_patches + num_flash_patches:batch_size, :]

            mprage_idx = np.random.randint(0, self.theta_mprage.shape[0])
            curr_mprage_theta = self.theta_mprage[mprage_idx, :]

            xmprage = self.apply_mprage(curr_mprage_theta, xmprage_pd, xmprage_t1, xmprage_t2)
            xmprage_feat = np.stack((xmprage,
                                     x_talx[num_orig_patches + num_flash_patches: batch_size, :],
                                     x_taly[num_orig_patches + num_flash_patches: batch_size, :],
                                     x_talz[num_orig_patches + num_flash_patches: batch_size, :],), axis=-2)

            xmprage_feat = xmprage_feat.reshape(xmprage_feat.shape[:-1])
            # print(xmprage_feat.shape)
            #
            # print(time.time() - t4)

            # stack xfeat, xflashfeat, xmpragefeat
            xallfeat = np.concatenate((xfeat, xflash_feat, xmprage_feat))
            x = xallfeat

            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_validation_generator_multichannel_nmr(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            # random.shuffle(rare_index_list)
            xfeat = np.stack((self.data_storage.root.val_srct1w[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctalx[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctaly[index_list[0:batch_size], :],
                              self.data_storage.root.val_srctalz[index_list[0:batch_size], :]), axis=-2)
            xfeat = xfeat.reshape(xfeat.shape[:-1])
            y = self.data_storage.root.val_trgseg[index_list[0:batch_size], :]
            x = xfeat
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_training_generator_multichannel_nmr_finetune(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        # num_orig_patches = batch_size//2
        num_flash_patches = batch_size
        # num_mprage_patches = batch_size - (num_orig_patches + num_flash_patches)
        while True:
            x_list = list()
            y_list = list()
            # t= time.time()
            random.shuffle(index_list)
            # print(time.time() - t)

            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalx[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctaly[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalz[index_list[0:num_orig_patches],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t2 = time.time()
            # x_t1w =  self.data_storage.root.srct1w[index_list[0:batch_size],:]
            x_talx = self.data_storage.root.srctalx[index_list[0:batch_size], :]
            x_taly = self.data_storage.root.srctaly[index_list[0:batch_size], :]
            x_talz = self.data_storage.root.srctalz[index_list[0:batch_size], :]

            x_pd = self.data_storage.root.srcnmrpd[index_list[0:batch_size], :]
            x_t1 = self.data_storage.root.srcnmrt1[
                   index_list[0:batch_size], :]
            x_t2 = self.data_storage.root.srcnmrt2[
                   index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]

            # print(time.time() - t2)
            # xfeat = np.stack((x_t1w[0:num_orig_patches,:],
            #                   x_talx[0:num_orig_patches,:],
            #                   x_taly[0:num_orig_patches, :],
            #                   x_talz[0:num_orig_patches, :]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t3 = time.time()
            xflash_pd = x_pd
            xflash_t1 = x_t1
            xflash_t2 = x_t2

            flash_idx = np.random.randint(0, self.theta_flash.shape[0])
            curr_theta_flash = self.theta_flash[flash_idx, :]

            xflash = self.apply_flash(curr_theta_flash, xflash_pd, xflash_t1, xflash_t2)
            xflash_feat = np.stack((xflash,
                                    x_talx,
                                    x_taly,
                                    x_talz,
                                    ), axis=-2)
            xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # print(time.time() - t3)
            # xflash_feat = np.stack((xflash,
            #                   self.data_storage.root.srctalx[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctaly[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctalz[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :]), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # # t4 = time.time()
            # xmprage_pd = x_pd[num_orig_patches + num_flash_patches:batch_size,:]
            # xmprage_t1 = x_t1[num_orig_patches + num_flash_patches:batch_size,:]
            # xmprage_t2 = x_t2[num_orig_patches + num_flash_patches:batch_size, :]

            # mprage_idx = np.random.randint(0, self.theta_mprage.shape[0])
            # curr_mprage_theta = self.theta_mprage[mprage_idx,:]
            #
            # xmprage = self.apply_mprage(curr_mprage_theta, xmprage_pd, xmprage_t1, xmprage_t2)
            #
            #
            #
            # xmprage_feat = np.stack((xmprage,
            #                          x_talx[num_orig_patches + num_flash_patches : batch_size,:],
            #                          x_taly[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_talz[num_orig_patches + num_flash_patches: batch_size, :],), axis=-2)
            #
            # xmprage_feat = xmprage_feat.reshape(xmprage_feat.shape[:-1])
            # print(xmprage_feat.shape)
            #
            # print(time.time() - t4)

            # stack xfeat, xflashfeat, xmpragefeat
            # xallfeat = np.concatenate((xfeat, xflash_feat, xmprage_feat))

            x = xflash_feat

            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_validation_generator_multichannel_nmr_finetune(self, batch_size):
        index_list = list(range(self.data_storage.root.val_srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        # num_orig_patches = batch_size//2
        num_flash_patches = batch_size
        # num_mprage_patches = batch_size - (num_orig_patches + num_flash_patches)
        while True:
            x_list = list()
            y_list = list()
            # t= time.time()
            random.shuffle(index_list)
            # print(time.time() - t)

            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalx[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctaly[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalz[index_list[0:num_orig_patches],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t2 = time.time()
            # x_t1w =  self.data_storage.root.srct1w[index_list[0:batch_size],:]
            x_talx = self.data_storage.root.srctalx[index_list[0:batch_size], :]
            x_taly = self.data_storage.root.srctaly[index_list[0:batch_size], :]
            x_talz = self.data_storage.root.srctalz[index_list[0:batch_size], :]

            x_pd = self.data_storage.root.srcnmrpd[index_list[0:batch_size], :]
            x_t1 = self.data_storage.root.srcnmrt1[
                   index_list[0:batch_size], :]
            x_t2 = self.data_storage.root.srcnmrt2[
                   index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]

            # print(time.time() - t2)
            # xfeat = np.stack((x_t1w[0:num_orig_patches,:],
            #                   x_talx[0:num_orig_patches,:],
            #                   x_taly[0:num_orig_patches, :],
            #                   x_talz[0:num_orig_patches, :]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t3 = time.time()
            xflash_pd = x_pd
            xflash_t1 = x_t1
            xflash_t2 = x_t2

            flash_idx = np.random.randint(0, self.theta_flash.shape[0])
            curr_theta_flash = self.theta_flash[flash_idx, :]

            xflash = self.apply_flash(curr_theta_flash, xflash_pd, xflash_t1, xflash_t2)
            xflash_feat = np.stack((xflash,
                                    x_talx,
                                    x_taly,
                                    x_talz,
                                    ), axis=-2)
            xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # print(time.time() - t3)
            # xflash_feat = np.stack((xflash,
            #                   self.data_storage.root.srctalx[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctaly[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctalz[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :]), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # # t4 = time.time()
            # xmprage_pd = x_pd[num_orig_patches + num_flash_patches:batch_size,:]
            # xmprage_t1 = x_t1[num_orig_patches + num_flash_patches:batch_size,:]
            # xmprage_t2 = x_t2[num_orig_patches + num_flash_patches:batch_size, :]

            # mprage_idx = np.random.randint(0, self.theta_mprage.shape[0])
            # curr_mprage_theta = self.theta_mprage[mprage_idx,:]
            #
            # xmprage = self.apply_mprage(curr_mprage_theta, xmprage_pd, xmprage_t1, xmprage_t2)
            #
            #
            #
            # xmprage_feat = np.stack((xmprage,
            #                          x_talx[num_orig_patches + num_flash_patches : batch_size,:],
            #                          x_taly[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_talz[num_orig_patches + num_flash_patches: batch_size, :],), axis=-2)
            #
            # xmprage_feat = xmprage_feat.reshape(xmprage_feat.shape[:-1])
            # print(xmprage_feat.shape)
            #
            # print(time.time() - t4)

            # stack xfeat, xflashfeat, xmpragefeat
            # xallfeat = np.concatenate((xfeat, xflash_feat, xmprage_feat))

            x = xflash_feat

            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_training_generator_augment(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_rare_patches = batch_size // 4
        num_remaining_patches = batch_size - num_rare_patches
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            random.shuffle(rare_index_list)
            for index in rare_index_list[:num_remaining_patches]:
                x_list.append(self.data_storage.root.src[index])
                y_list.append(self.data_storage.root.trg[index])

            for index in rare_index_list[:num_rare_patches]:
                x_list.append(self.data_storage.root.src_rare[index])
                y_list.append(self.data_storage.root.trg_rare[index])

            # probabilistically augment (rotate a batch
            # x, y = self.convert_data(x_list, y_list, n_labels=self.n_labels, labels=self.labels)
            x = np.asarray(x_list)
            y = np.asarray(y_list)

            # num_rotated_patches = int(batch_size * self.prob_rot_augmentation)
            # rot_idxs = np.random.choice(batch_size, size=(num_rotated_patches,), replace=False)
            # all_axes = ((0,1), (1,2), (2,0))
            # # print(rot_idxs)
            # for r_idx in rot_idxs:
            #     # choose angle between uniformly randomly between -30 and 30
            #     r_angle = 60*np.random.rand() - 30
            #     ax_idx = np.random.choice(3, size=(1,), replace=False)
            #     axes = all_axes[ax_idx[0]]
            #
            #
            #     rx = scipy.ndimage.interpolation.rotate(x[r_idx,:,:,:,0], angle=r_angle, axes=axes, reshape=False, mode='nearest')
            #     rx[rx<0] = 0
            #     ry = scipy.ndimage.interpolation.rotate(x[r_idx,:,:,:,0], angle=r_angle, axes=axes, order=0, reshape=False, mode='nearest')
            #     ry[ry<0] = 0
            #
            #     x[r_idx,:,:,:,0] = rx
            #     y[r_idx, :, :, :, 0] = ry

            # probabilistically select half of the patches in x and add a intensity augmentation
            num_int_patches = int(batch_size * self.prob_int_augmentation)
            int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # print(int_idxs)
            for i_idx in int_idxs:
                intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
                x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_training_generator_singlechannel_nmr(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_orig_patches = batch_size // 2
        num_flash_patches = batch_size // 4
        num_mprage_patches = batch_size - (num_orig_patches + num_flash_patches)

        while True:
            x_list = list()
            y_list = list()
            # t= time.time()
            random.shuffle(index_list)
            # print(time.time() - t)

            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalx[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctaly[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalz[index_list[0:num_orig_patches],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t2 = time.time()
            x_t1w = self.data_storage.root.srct1w[index_list[0:batch_size], :]
            # x_talx = self.data_storage.root.srctalx[index_list[0:batch_size], :]
            # x_taly = self.data_storage.root.srctaly[index_list[0:batch_size], :]
            # x_talz = self.data_storage.root.srctalz[index_list[0:batch_size], :]

            x_pd = self.data_storage.root.srcnmrpd[index_list[0:batch_size], :]
            x_t1 = self.data_storage.root.srcnmrt1[
                   index_list[0:batch_size], :]
            x_t2 = self.data_storage.root.srcnmrt2[
                   index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]

            # print(time.time() - t2)

            # xfeat = np.stack((x_t1w[0:num_orig_patches, :],
            #                   x_talx[0:num_orig_patches, :],
            #                   x_taly[0:num_orig_patches, :],
            #                   x_talz[0:num_orig_patches, :]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t3 = time.time()
            xflash_pd = x_pd[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t1 = x_t1[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t2 = x_t2[num_orig_patches:num_orig_patches + num_flash_patches, :]

            flash_idx = np.random.randint(0, self.theta_flash.shape[0])
            curr_theta_flash = self.theta_flash[flash_idx, :]

            xflash = self.apply_flash(curr_theta_flash, xflash_pd, xflash_t1, xflash_t2)
            # xflash_feat = np.stack((xflash,
            #                         x_talx[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         x_taly[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         x_talz[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         ), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # print(time.time() - t3)

            # xflash_feat = np.stack((xflash,
            #                   self.data_storage.root.srctalx[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctaly[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctalz[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :]), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # t4 = time.time()
            xmprage_pd = x_pd[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t1 = x_t1[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t2 = x_t2[num_orig_patches + num_flash_patches:batch_size, :]

            mprage_idx = np.random.randint(0, self.theta_mprage.shape[0])
            curr_mprage_theta = self.theta_mprage[mprage_idx, :]

            xmprage = self.apply_mprage(curr_mprage_theta, xmprage_pd, xmprage_t1, xmprage_t2)

            # xmprage_feat = np.stack((xmprage,
            #                          x_talx[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_taly[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_talz[num_orig_patches + num_flash_patches: batch_size, :],), axis=-2)

            # xmprage_feat = xmprage_feat.reshape(xmprage_feat.shape[:-1])
            # print(xmprage_feat.shape)
            #
            # print(time.time() - t4)

            # stack xfeat, xflashfeat, xmpragefeat
            xallfeat = np.concatenate((x_t1w, xflash, xmprage))
            # xallfeat = np.stack((x, xflash, xmprage), axis=-2)
            # xallfeat.reshape(xallfeat.shape[:-1])

            x = xallfeat

            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def seg_validation_generator_singlechannel_nmr(self, batch_size):
        index_list = list(range(self.data_storage.root.srct1w.shape[0]))
        # rare_index_list = list(range(self.data_storage.root.src_rare.shape[0]))
        num_orig_patches = batch_size // 2
        num_flash_patches = batch_size // 4
        num_mprage_patches = batch_size - (num_orig_patches + num_flash_patches)
        while True:
            x_list = list()
            y_list = list()
            # t= time.time()
            random.shuffle(index_list)
            # print(time.time() - t)

            # xfeat = np.stack((self.data_storage.root.srct1w[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalx[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctaly[index_list[0:num_orig_patches],:],
            #                   self.data_storage.root.srctalz[index_list[0:num_orig_patches],:]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])

            # t2 = time.time()
            x_t1w = self.data_storage.root.srct1w[index_list[0:batch_size], :]
            # x_talx = self.data_storage.root.srctalx[index_list[0:batch_size], :]
            # x_taly = self.data_storage.root.srctaly[index_list[0:batch_size], :]
            # x_talz = self.data_storage.root.srctalz[index_list[0:batch_size], :]

            x_pd = self.data_storage.root.srcnmrpd[index_list[0:batch_size], :]
            x_t1 = self.data_storage.root.srcnmrt1[
                   index_list[0:batch_size], :]
            x_t2 = self.data_storage.root.srcnmrt2[
                   index_list[0:batch_size], :]
            y = self.data_storage.root.trgseg[index_list[0:batch_size], :]

            # print(time.time() - t2)

            # xfeat = np.stack((x_t1w[0:num_orig_patches, :],
            #                   x_talx[0:num_orig_patches, :],
            #                   x_taly[0:num_orig_patches, :],
            #                   x_talz[0:num_orig_patches, :]), axis=-2)
            # xfeat = xfeat.reshape(xfeat.shape[:-1])
            x_t1w = x_t1w[0:num_orig_patches, :]
            # t3 = time.time()
            xflash_pd = x_pd[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t1 = x_t1[num_orig_patches:num_orig_patches + num_flash_patches, :]
            xflash_t2 = x_t2[num_orig_patches:num_orig_patches + num_flash_patches, :]

            flash_idx = np.random.randint(0, self.theta_flash.shape[0])
            curr_theta_flash = self.theta_flash[flash_idx, :]

            xflash = self.apply_flash(curr_theta_flash, xflash_pd, xflash_t1, xflash_t2)
            # xflash_feat = np.stack((xflash,
            #                         x_talx[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         x_taly[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         x_talz[num_orig_patches:num_orig_patches + num_flash_patches],
            #                         ), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # print(time.time() - t3)
            # xflash_feat = np.stack((xflash,
            #                   self.data_storage.root.srctalx[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctaly[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :],
            #                   self.data_storage.root.srctalz[index_list[num_orig_patches:num_orig_patches + num_flash_patches], :]), axis=-2)
            # xflash_feat = xflash_feat.reshape(xflash_feat.shape[:-1])

            # t4 = time.time()
            xmprage_pd = x_pd[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t1 = x_t1[num_orig_patches + num_flash_patches:batch_size, :]
            xmprage_t2 = x_t2[num_orig_patches + num_flash_patches:batch_size, :]

            mprage_idx = np.random.randint(0, self.theta_mprage.shape[0])
            curr_mprage_theta = self.theta_mprage[mprage_idx, :]

            xmprage = self.apply_mprage(curr_mprage_theta, xmprage_pd, xmprage_t1, xmprage_t2)

            # xmprage_feat = np.stack((xmprage,
            #                          x_talx[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_taly[num_orig_patches + num_flash_patches: batch_size, :],
            #                          x_talz[num_orig_patches + num_flash_patches: batch_size, :],), axis=-2)

            # xmprage_feat = xmprage_feat.reshape(xmprage_feat.shape[:-1])
            # print(xmprage_feat.shape)
            #
            # print(time.time() - t4)

            # stack xfeat, xflashfeat, xmpragefeat
            xallfeat = np.concatenate((x_t1w, xflash, xmprage))
            print(xallfeat.shape)
            print(y.shape)
            # xallfeat = np.stack((x, xflash, xmprage), axis=-2)
            # xallfeat.reshape(xallfeat.shape[:-1])

            x = xallfeat

            # num_int_patches = int(batch_size * self.prob_int_augmentation)
            # int_idxs = np.random.choice(batch_size, size=(num_int_patches,), replace=False)
            # # print(int_idxs)
            # for i_idx in int_idxs:
            #     intx = self.perturb_intensity(x[i_idx, :, :, :, 0])
            #     x[i_idx, :, :, :, 0] = intx

            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def perturb_intensity(self, patch):
        # randomly choose an integer between 0 and 51
        hist_idx = np.random.randint(0, self.augment_cubic_transforms.shape[0] - 1, 1)[0]
        cubic_transform = self.augment_cubic_transforms[hist_idx, :]
        poly_cubic = np.poly1d(cubic_transform)
        aug_patch = poly_cubic(patch)
        # interp_ref_values = self.interp_ref_values[:,hist_idx]
        # bins_in = bins_in_z[1:]
        # aug_patch = np.zeros(patch.shape)
        # for i in range(1, len(bins_in)):
        #     aug_patch[(patch > bins_in_z[i - 1]) & (patch <= bins_in_z[i])] = interp_ref_values[i - 1]
        return aug_patch

    def seg_validation_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src_validation.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src_validation[index])
                y_list.append(self.data_storage.root.trg_validation[index])
            x = np.asarray(x_list)
            y = np.asarray(y_list)
            x, y = self.convert_data(x, y, n_labels=self.n_labels, labels=self.labels)
            yield x, y

    def training_label_generator(self, batch_size):
        index_list = list(range(self.data_storage.root.src.shape[0]))
        while True:
            x_list = list()
            y_list = list()
            random.shuffle(index_list)
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
            random.shuffle(index_list)
            for index in index_list[:batch_size]:
                x_list.append(self.data_storage.root.src_validation[index])
                y_list.append(self.data_storage.root.trg_validation[index])
            x_list = np.asarray(x_list)
            y_list = np.asarray(y_list)
            yield x_list, y_list
