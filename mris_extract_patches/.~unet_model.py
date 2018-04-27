from __future__ import print_function

from keras.models import Model
from keras.layers import Input, concatenate, \
    Conv2D, MaxPooling2D, Conv2DTranspose, Conv3D,MaxPooling3D, Activation,\
    Deconvolution3D,UpSampling3D,UpSampling2D, BatchNormalization, ZeroPadding3D, Cropping3D, MaxPool3D
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint
from keras import backend as K
from keras.losses import  mean_squared_error, mean_absolute_percentage_error
from keras.models import Sequential
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten
import keras
from keras.utils.training_utils import  multi_gpu_model
import numpy as np
import os
import tensorflow as tf
import SBRegCF
from SBRegCF import SBRegCF


from keras.utils.generic_utils import get_custom_objects

K.set_image_data_format('channels_last')


def dice_coef_loss2(y_true, y_pred):
    area_reg = 0.1
    y_true /= K.sum(y_true, axis=-1, keepdims=True)
    y_true = K.clip(y_true, K.epsilon(), 1)

    # make sure pred is a probability
    y_pred /= K.sum(y_pred, axis=-1, keepdims=True)
    y_pred = K.clip(y_pred, K.epsilon(), 1)
    y_pred_op = y_pred
    y_true_op = y_true


    # compute dice for each entry in batch.
    # dice will now be [batch_size, nb_labels]
    sum_dim = 1
    top = 2 * K.sum(y_true_op * y_pred_op, sum_dim)
    bottom = K.sum(K.square(y_true_op), sum_dim) + K.sum(K.square(y_pred_op), sum_dim)
    # make sure we have no 0s on the bottom. K.epsilon()
    bottom = K.maximum(bottom, area_reg)
    dice_metric = top / bottom
    dice_loss = 1 - dice_metric
    mean_dice_loss = K.mean(dice_loss)

    return mean_dice_loss




def dice_coef(y_true, y_pred, smooth=1.):
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)


def dice_coef_loss(y_true, y_pred):
    return -dice_coef(y_true, y_pred)


def pure_grad_loss(y_true, y_pred):
    print('grad loss')
    p0 = y_true.shape[1]
    p1 = y_true.shape[2]
    p2 = y_true.shape[3]
    lambda1 = 10
    sobel_z = np.zeros((3,3,3))
    sobel_z[:, :, 0] = [[-1.0, -3.0, -1.0], [-3.0, -6.0, -3.0], [-1.0, -3.0, -1.0]]
    sobel_z[:, :, 1] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    sobel_z[:, :, 2] = -sobel_z[:, :, 0]
    sobel_z = sobel_z/22
    sobel_x = np.swapaxes(sobel_z, 0, 2)
    sobel_y = np.swapaxes(sobel_z, 1, 2)
    sobel_z = np.reshape(sobel_z, (3, 3, 3, 1, 1))
    sobel_y = np.reshape(sobel_y, (3, 3, 3, 1, 1))
    sobel_x = np.reshape(sobel_x, (3, 3, 3, 1, 1))




    tf_sobel_z = tf.convert_to_tensor(sobel_z, dtype='float32', name='sobel_z')
    tf_sobel_y = tf.convert_to_tensor(sobel_y, dtype='float32', name='sobel_y')
    tf_sobel_x = tf.convert_to_tensor(sobel_x, dtype='float32', name='sobel_x')

    y_true_zgrad = tf.nn.conv3d(y_true, tf_sobel_z, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_true_ygrad = tf.nn.conv3d(y_true, tf_sobel_y, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_true_xgrad = tf.nn.conv3d(y_true, tf_sobel_x, strides=[1, 1, 1, 1, 1], padding='SAME')

    y_pred_zgrad = tf.nn.conv3d(y_pred, tf_sobel_z, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_pred_ygrad = tf.nn.conv3d(y_pred, tf_sobel_y, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_pred_xgrad = tf.nn.conv3d(y_pred, tf_sobel_x, strides=[1, 1, 1, 1, 1], padding='SAME')



    return lambda1 * mean_squared_error(y_true_zgrad, y_pred_zgrad) +\
           lambda1 * mean_squared_error(y_true_ygrad, y_pred_ygrad) +\
           lambda1 * mean_squared_error(y_true_xgrad, y_pred_xgrad)


def grad_loss(y_true, y_pred):
    print('grad loss')
    p0 = y_true.shape[1]
    p1 = y_true.shape[2]
    p2 = y_true.shape[3]
    lambda1 = 10000
    sobel_z = np.zeros((3,3,3))
    sobel_z[:, :, 0] = [[-1.0, -3.0, -1.0], [-3.0, -6.0, -3.0], [-1.0, -3.0, -1.0]]
    sobel_z[:, :, 1] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    sobel_z[:, :, 2] = -sobel_z[:, :, 0]
    sobel_z = sobel_z/22
    sobel_x = np.swapaxes(sobel_z, 0, 2)
    sobel_y = np.swapaxes(sobel_z, 1, 2)
    sobel_z = np.reshape(sobel_z, (3, 3, 3, 1, 1))
    sobel_y = np.reshape(sobel_y, (3, 3, 3, 1, 1))
    sobel_x = np.reshape(sobel_x, (3, 3, 3, 1, 1))




    tf_sobel_z = tf.convert_to_tensor(sobel_z, dtype='float32', name='sobel_z')
    tf_sobel_y = tf.convert_to_tensor(sobel_y, dtype='float32', name='sobel_y')
    tf_sobel_x = tf.convert_to_tensor(sobel_x, dtype='float32', name='sobel_x')

    y_true_zgrad = tf.nn.conv3d(y_true, tf_sobel_z, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_true_ygrad = tf.nn.conv3d(y_true, tf_sobel_y, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_true_xgrad = tf.nn.conv3d(y_true, tf_sobel_x, strides=[1, 1, 1, 1, 1], padding='SAME')

    y_pred_zgrad = tf.nn.conv3d(y_pred, tf_sobel_z, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_pred_ygrad = tf.nn.conv3d(y_pred, tf_sobel_y, strides=[1, 1, 1, 1, 1], padding='SAME')
    y_pred_xgrad = tf.nn.conv3d(y_pred, tf_sobel_x, strides=[1, 1, 1, 1, 1], padding='SAME')



    return mean_squared_error(y_true, y_pred) + \
           lambda1 * mean_squared_error(y_true_zgrad, y_pred_zgrad) +\
           lambda1 * mean_squared_error(y_true_ygrad, y_pred_ygrad) +\
           lambda1 * mean_squared_error(y_true_xgrad, y_pred_xgrad)


def unet_model_3d(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=0,
                  loss='mean_absolute_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True, num_gpus=1):
    """
    Builds the 3D UNet Keras model.
    :param input_shape: Shape of the input data (x_size, y_size, z_size).
    :param downsize_filters_factor: Factor to which to reduce the number of filters. Making this value larger will
    reduce the amount of memory the model will need during training.
    :param pool_size: Pool size for the max pooling operations.
    :param n_labels: Number of binary labels that the model is learning.
    :param initial_learning_rate: Initial learning rate for the model. This will be decayed during training.
    :param deconvolution: If set to True, will use transpose convolution(deconvolution) instead of upsamping. This
    increases the amount memory required during training.
    :return: Untrained 3D UNet Model
    """
    # channels last, make feature shape from (32,32,32) - > (32,32,32,1)

    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []
    print('unet depth is ')
    print(unet_depth)
    for i in range(unet_depth):

        prev = input_img if i == 0 else pools[i-1]
        print(int(num_filters*(2**i)/downsize_filters_factor))
        conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_1_%d' % (i)))(prev)
        conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
        conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_2_%d' % (i)))(conv)
        conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
        if i < (unet_depth - 1):
            pools.append(MaxPooling3D(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))

        convs.append(conv)

    for i in range(unet_depth - 1):
        index = i + unet_depth - 1
        level = unet_depth - (i + 2)
        up = concatenate([UpSampling3D(size=pool_size,  name=('upsampling_U_%d' % (level+1)))(convs[index]),
                          convs[level]], axis=-1,  name=('concat_%d' % (level)))
        conv = Conv3D(num_filters * (2 ** level), (3, 3, 3), padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv3D_U_1_%d' % (level))
                      )(up)
        conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
        conv = Conv3D(num_filters * (2 ** level), (3, 3, 3), padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv3D_U_2_%d' % (level)))(conv)
        convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))

#    conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
#    conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu",
#                  kernel_initializer="he_normal")(conv)
#    conv = BatchNormalization()(conv)
#    center_input = Cropping3D(cropping=(0, 0, 0))(input_img)

    inputs.append(input_img)
#    centered_inputs.append(center_input)
    print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    print(loss)
    print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        conv = Conv3D(1, (1,1,1), activation='relu',  name='final_conv_3d')(up)


        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=inputs, outputs=conv)
                parallel_model = multi_gpu_model(model, gpus=num_gpus)
                if loss == 'grad_loss':
                    print(loss)
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                else:
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)

                return model, parallel_model
        else:
            model = Model(inputs=inputs, outputs=conv)
            if loss == 'grad_loss':
                print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        print('segmentation network')
        if n_labels > 1:
            conv = Conv3D(n_labels, (1, 1, 1), activation='softmax', name='final_conv_3d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model
        else:
            conv = Conv3D(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model


def unet_model_2d_old(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2), n_labels=0,
                  loss='mean_absolute_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True, 
                  num_gpus=1, num_outputs=1):
    """
    Builds the 2D UNet Keras model.
    :param input_shape: Shape of the input data (x_size, y_size).
    :param downsize_filters_factor: Factor to which to reduce the number of filters. Making this value larger will
    reduce the amount of memory the model will need during training.
    :param pool_size: Pool size for the max pooling operations.
    :param n_labels: Number of binary labels that the model is learning.
    :param initial_learning_rate: Initial learning rate for the model. This will be decayed during training.
    :param deconvolution: If set to True, will use transpose convolution(deconvolution) instead of upsamping. This
    increases the amount memory required during training.
    :return: Untrained 2D UNet Model
    """
    # channels last, make feature shape from (32,32,32) - > (32,32,32,1)

    ConvL = Conv2D
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input' )
    convs = []
    pools = []
    inputs = []
    out_filter_shape = (1,1)
    centered_inputs = []
    endpoints = []
    print('unet depth is ')
    print(unet_depth)
    for i in range(unet_depth):

        prev = input_img if i == 0 else pools[i-1]
        print(int(num_filters*(2**i)/downsize_filters_factor))
        conv = Conv2D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv2D_D_1_%d' % (i)))(prev)
        conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
        conv = Conv2D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv2D_D_2_%d' % (i)))(conv)
        conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
        if i < (unet_depth - 1):
            pools.append(MaxPooling2D(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))

        convs.append(conv)

    for i in range(unet_depth - 1):
        index = i + unet_depth - 1
        level = unet_depth - (i + 2)
        print('index is ')
        print(index)
        up = concatenate([UpSampling2D(size=pool_size,  name=('upsampling_U_%d' % (level+1)))(convs[index]),
                          convs[level]], axis=-1,  name=('concat_%d' % (level)))
        conv = Conv2D(num_filters * (2 ** level), (3, 3), padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv2D_U_1_%d' % (level))
                      )(up)
        conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
        conv = Conv2D(num_filters * (2 ** level), (3, 3), padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv2D_U_2_%d' % (level)))(conv)
        convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))

#    conv = ZeroPadding2D(padding=(1, 1, 1))(convs[-1])
#    conv = Conv2D(num_filters * 2, (3, 3), padding="valid", activation="relu",
#                  kernel_initializer="he_normal")(conv)
#    conv = BatchNormalization()(conv)
#    center_input = Cropping2D(cropping=(0, 0, 0))(input_img)

    inputs.append(input_img)
#    centered_inputs.append(center_input)
    print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    print(loss)
    print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        conv = ConvL(num_outputs, out_filter_shape, activation='relu',  name='final_conv_2d')(up)
#        conv = Conv2D(1, (1,1), activation='relu',  name='final_conv_2d')(up)


        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=inputs, outputs=conv)
                parallel_model = multi_gpu_model(model, gpus=num_gpus)
                if loss == 'grad_loss':
                    print(loss)
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                else:
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)

                return model, parallel_model
        else:
            model = Model(inputs=inputs, outputs=conv)
            if loss == 'grad_loss':
                print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        print('segmentation network')
        if n_labels > 1:
            conv = Conv2D(n_labels, (1, 1, 1), activation='softmax', name='final_conv_2d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model
        else:
            conv = Conv2D(1, (1, 1, 1), activation='sigmoid', name='final_conv_2d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model


def class_net(feature_shape, dim, unet_num_filters, n_labels, initial_learning_rate, loss='binary_cross_entropy'):
    if dim == 2:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        filter_shape = (3,3)
        pool_shape = (2,2)
    elif dim ==3:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        filter_shape = (3,3,3)
        pool_shape = (2,2,2)

    input_shape_list = list(feature_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    model = Sequential()
    model.add(ConvL(unet_num_filters, filter_shape, padding='same', activation='relu', input_shape=input_shape_append))
    model.add(ConvL(unet_num_filters, filter_shape, padding='valid', activation='relu'))
    model.add(MaxPoolingL(pool_size=pool_shape))
    model.add(Dropout(0.25))

    model.add(ConvL(2*unet_num_filters, filter_shape, padding='same', activation='relu'))
    model.add(ConvL(2*unet_num_filters, filter_shape, activation='relu'))
    model.add(MaxPoolingL(pool_size=pool_shape))
    model.add(Dropout(0.25))

    model.add(ConvL((2**2)*unet_num_filters, filter_shape, padding='same', activation='relu'))
    model.add(ConvL((2**2)*unet_num_filters, filter_shape, activation='relu'))
    model.add(MaxPoolingL(pool_size=pool_shape))
    model.add(Dropout(0.25))

    model.add(ConvL((2**3)*unet_num_filters, filter_shape, padding='same', activation='relu'))
    model.add(ConvL((2**3)*unet_num_filters, filter_shape, activation='relu'))
    model.add(MaxPoolingL(pool_size=pool_shape))
    model.add(Dropout(0.25))

    model.add(ConvL((2**4)*unet_num_filters, filter_shape, padding='same', activation='relu'))
    model.add(ConvL((2**4)*unet_num_filters, filter_shape, activation='relu'))
    model.add(MaxPoolingL(pool_size=pool_shape))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(512, activation='relu'))
    model.add(Dropout(0.5))
    if n_labels == 2:
        model.add(Dense(1, activation='sigmoid'))
    else:
        model.add(Dense(n_labels, activation='softmax'))
    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss, metrics=['accuracy'])
    return model


def atrous_net(input_shape, num_filters,initial_learning_rate=0.00001, loss='mean_absolute_error'):

    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []

    x = Conv3D(num_filters, (3, 3, 3), activation='relu', padding='same', dilation_rate=1)(input_img)
    x = Conv3D(num_filters, (3, 3, 3), activation='relu', padding='same', dilation_rate=1)(x)
    x = Conv3D(num_filters, (3, 3, 3), activation='relu', padding='same', dilation_rate=2)(x)
    x = Conv3D(num_filters, (3, 3, 3), activation='relu', padding='same', dilation_rate=4)(x)
    x = Conv3D(1, (1, 1, 1), activation='relu', name='final_conv_3d')(x)

    model = Model(inputs=input_img, outputs=x)
    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
    return model


def noise_net(input_shape, num_filters,initial_learning_rate=0.00001, loss='mean_absolute_error'):
    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []

    x = Conv3D(num_filters, (7, 7, 7), activation='relu', padding='same', dilation_rate=1)(input_img)
    x = Conv3D(num_filters/2, (5,5, 5), activation='relu', padding='same', dilation_rate=1)(x)
    x = Conv3D(1, (1, 1, 1), activation='relu', name='final_conv_3d')(x)
    model = Model(inputs=input_img, outputs=x)
    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
    return model


def compute_level_output_shape(filters, depth, pool_size, image_shape):
    """
    Each level has a particular output shape based on the number of filters used in that level and the depth or number
    of max pooling operations that have been done on the data at that point.
    :param image_shape: shape of the 3d image.
    :param pool_size: the pool_size parameter used in the max pooling operation.
    :param filters: Number of filters used by the last node in a given level.
    :param depth: The number of levels down in the U-shaped model a given node is.
    :return: 5D vector of the shape of the output node
    """
    if depth != 0:
        output_image_shape = np.divide(image_shape, np.multiply(pool_size, depth)).tolist()
    else:
        output_image_shape = image_shape
    return tuple([None, filters] + [int(x) for x in output_image_shape])




def get_upconv(depth, nb_filters, pool_size, image_shape, kernel_size=(2, 2, 2), strides=(2, 2, 2),
               deconvolution=False):
    if deconvolution:
        return Deconvolution3D(filters=nb_filters, kernel_size=kernel_size,
                               output_shape=compute_level_output_shape(filters=nb_filters, depth=depth,
                                                                       pool_size=pool_size, image_shape=image_shape),
                               strides=strides, input_shape=compute_level_output_shape(filters=nb_filters,
                                                                                       depth=depth+1,
                                                                                       pool_size=pool_size,
                                                                                       image_shape=image_shape))
    else:
        return UpSampling3D(size=pool_size)

if __name__=="__main__":
    os.environ['LD_LIBRARY_PATH'] = '/usr/pubsw/packages/CUDA/lib64 '
    # model = class_net(feature_shape=(256,256), dim=2, unet_num_filters=32, unet_depth=4, n_labels=2,
    #                   initial_learning_rate=0.001, loss='categorical_crossentropy')
    # model = atrous_net(input_shape=(16, 16, 16), num_filters=32)
    # model = unet_model_3d(input_shape = (16,16,16), num_filters=32, unet_depth=2,
    #                       downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=10,
    #                       loss='dice_coef_loss', initial_learning_rate=0.00001, deconvolution=False)


def unet_model_2d(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2), n_labels=0,
                  loss='mean_squared_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True,
                  num_gpus=1, num_outputs=1):
    """
    Builds the 3D UNet Keras model.
    :param input_shape: Shape of the input data (x_size, y_size, z_size).
    :param downsize_filters_factor: Factor to which to reduce the number of filters. Making this value larger will
    reduce the amount of memory the model will need during training.
    :param pool_size: Pool size for the max pooling operations.
    :param n_labels: Number of binary labels that the model is learning.
    :param initial_learning_rate: Initial learning rate for the model. This will be decayed during training.
    :param deconvolution: If set to True, will use transpose convolution(deconvolution) instead of upsamping. This
    increases the amount memory required during training.
    :return: Untrained 3D UNet Model
    """
    # channels last, make feature shape from (32,32,32) - > (32,32,32,1)

    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    dim = len(input_shape)
    if dim == 3:
        ConvL =  Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2,2)
        UpSamplingL = UpSampling2D
        filter_shape = (5,5)
        out_filter_shape = (1,1)

    elif dim==4:
        ConvL = Conv3D
        MaxPoolingL= MaxPooling3D
        pool_size = (2,2,2)
        UpSamplingL = UpSampling3D
        filter_shape = (3,3,3)
        out_filter_shape = (1,1,1)

    print('out filter shape is ' + str(out_filter_shape))
    # input_shape_list = list(input_shape)
    # input_shape_list.append(1)
    # input_shape_append = tuple(input_shape_list)
    # print(input_shape_append)
    input_img = Input(shape=input_shape, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []

    print('unet depth is ')
    print(unet_depth)
    for i in range(unet_depth):

        prev = input_img if i == 0 else pools[i-1]
        print(int(num_filters*(2**i)/downsize_filters_factor))
        conv = ConvL(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_1_%d' % (i)))(prev)
        # conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
        conv = ConvL(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_2_%d' % (i)))(conv)
        conv = ConvL(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_1_%d' % (i)))(prev)
        if i < (unet_depth - 1):
            pools.append(MaxPoolingL(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))

        convs.append(conv)

    for i in range(unet_depth - 1):
        index = i + unet_depth - 1
        level = unet_depth - (i + 2)
        up = concatenate([UpSamplingL(size=pool_size,  name=('upsampling_U_%d' % (level+1)))(convs[index]),
                          convs[level]], axis=-1,  name=('concat_%d' % (level)))
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv3D_U_1_%d' % (level))
                      )(up)
        # conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                      kernel_initializer="he_normal",
                      name=('conv3D_U_2_%d' % (level)))(conv)
        # convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))
        convs.append(conv)

#    conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
#    conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu",
#                  kernel_initializer="he_normal")(conv)
#    conv = BatchNormalization()(conv)
#    center_input = Cropping3D(cropping=(0, 0, 0))(input_img)

    inputs.append(input_img)
#    centered_inputs.append(center_input)
    print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    print(loss)
    print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        conv = ConvL(num_outputs, out_filter_shape, activation='linear',  name='final_conv_3d')(up)
        conv = Flatten()(conv)


        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=inputs, outputs=conv)
                parallel_model = multi_gpu_model(model, gpus=num_gpus)
                if loss == 'grad_loss':
                    print(loss)
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                else:
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)

                return model, parallel_model
        else:
            model = Model(inputs=inputs, outputs=conv)
            if loss == 'grad_loss':
                print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            elif loss == 'warp_image_loss':
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=warp_image_loss)
            elif loss == 'warp_image_loss1':
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=warp_image_loss1)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        print('segmentation network')
        if n_labels > 1:
            conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model
        else:
            conv = Conv3D(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(up)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=inputs, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=inputs, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model


from keras import backend as K
import warp_image_utils
from warp_image_utils import warp_image
import tensorflow as tf
import nibabel as nib
from spatial_transformer_network import batch_warp2d, batch_mgrid, batch_displacement_warp2d
#from spatial_transformer_network import grid
import grid

def warp_jacobian(warp_tensor, pad=15):
    width = 512
    height = 512
    x0 = tf.linspace(0.0, width-1.0, width)
    y0 = tf.linspace(0.0, height-1.0, height)
    p0x, p0y = tf.meshgrid(x0, y0)
    
    p0xw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0x, (1, width, height, 1)), warp_tensor), (width, height))
    p0yw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0y, (1, width, height, 1)), warp_tensor), (width, height))
    
    p1xw = tf.manip.roll(p0xw, shift=[1,0], axis=[0,1]) ;
    p1yw = tf.manip.roll(p0yw, shift=[1,0], axis=[0,1]) ;
    
    p2xw = tf.manip.roll(p0xw, shift=[0,1], axis=[0,1]) ;
    p2yw = tf.manip.roll(p0yw, shift=[0,1], axis=[0,1]) ;
    
    de1x = p1xw - p0xw 
    de1y = p1yw - p0yw
    
    de2x = p2xw - p0xw 
    de2y = p2yw - p0yw

    jac = de2x*de1y - de2y*de1x
    pad0 = tf.cast(p0x > pad, tf.float32)
    pad1 = tf.cast(p0x < width-pad, tf.float32)
    pad2 = tf.cast(p0y > pad, tf.float32)
    pad3 = tf.cast(p0y < height-pad, tf.float32)
    pad = pad0*pad1*pad2*pad3
    return jac*pad

def compute_magnitude_loss(warp_tensor, threshold=512.0/512.0):
    dx = warp_tensor[0,0,:,:] 
    dy = warp_tensor[0,1,:,:] 
    zdx = dx - tf.reduce_mean(dx) ;
    zdy = dy - tf.reduce_mean(dy) ;
    norm = tf.sqrt(tf.add(dx ** 2, dy ** 2));
    mag_loss = tf.sqrt(tf.reduce_mean(((1-threshold)+norm)**2)) ;
#    mag_loss = tf.sqrt(tf.reduce_mean(norm**2)) ;
    return(mag_loss)

def warp_image_loss(y_true, y_pred):
    """
    applies warp (width x height x 2) to img (width x height) and returns the warped image
    first column of y_true is input image, and 2nd is desired output image
    """
    return compute_warp_image_loss(y_true, y_pred, image_weight=0, smooth_weight=1.0, jac_weight=1.,sb_weight=.01,mag_weight=10)

def correlation_loss(image1, image2):
    """
    computed the correlation of two images
    """
    z1 = image1 - tf.reduce_mean(image1)
    z2 = image2 - tf.reduce_mean(image2)
    corr = tf.reduce_mean(z1 * z2)
    return(corr)

def warp_image_loss1(y_true, y_pred):
    """
    applies warp (width x height x 2) to img (width x height) and returns the warped image
    first column of y_true is input image, and 2nd is desired output image
    """
    
    width = 512
    height = 512
    big_tensor = tf.reshape(y_true, (1, 4, width, height))
    images = big_tensor[:, 0:2, :,:]
    dx0 = tf.squeeze(big_tensor[:, 2,:,:])
    dy0 = tf.squeeze(big_tensor[:, 3,:,:])
    warp_tensor = tf.reshape(y_pred, (1, 2, width, height))/(512) ;
    zwarp_tensor = warp_tensor - tf.reshape(tf.stack([dx0, dy0]), (1, 2, width, height));
    mag_loss = compute_magnitude_loss(warp_tensor,30./width) ;
    loss = compute_warp_image_loss(images, y_pred, image_weight=1, smooth_weight=.0, jac_weight=0.0,sb_weight=0,mag_weight=0,mi_weight=0)
    real_image = tf.squeeze(tf.cast(images[0, 0, :, :], dtype=tf.float32))
    pred_image = tf.squeeze(tf.cast(images[0, 1, :, :], dtype=tf.float32))
    test = real_image * pred_image
    corr_loss = correlation_loss(real_image, pred_image) ;
    mag_loss *= 0.001*mag_loss
    return loss


def compute_warp_image_loss(y_true, y_pred, image_weight=1,smooth_weight=1,jac_weight=0,sb_weight=0, mag_weight = 0,mi_weight=0):
    """
    applies warp (width x height x 2) to img (width x height) and returns the warped image
    first column of y_true is input image, and 2nd is desired output image
    """
    width = 512
    height = 512
    smooth_iter = 100
    smooth_warp_iter = 0


    laplace_weight = 0

    images = tf.reshape(y_true, (1, 2*width, height))

    real_image = tf.cast(images[0, 0:width, :], dtype=tf.float32)
    pred_image = images[0, width:2 * width, :]
    x0 = tf.linspace(0.0, width-1.0, width)
    y0 = tf.linspace(0.0, height-1.0, height)
    p0x, p0y = tf.meshgrid(x0, y0)
    pad = 15;
    pad0 = tf.cast(p0x > pad, tf.float32)
    pad1 = tf.cast(p0x < width-pad, tf.float32)
    pad2 = tf.cast(p0y > pad, tf.float32)
    pad3 = tf.cast(p0y < height-pad, tf.float32)
    pad = pad0*pad1*pad2*pad3
    warp_tensor = tf.reshape(y_pred, (1, 2, width, height))/(512) ;

    dx = warp_tensor[0,0,:,:] 
    dy = warp_tensor[0,1,:,:] 
    warp_tensor = tf.reshape(tf.stack([dx, dy]), (1, 2, width, height))

    inverse_warp_tensor = warp_inverse(warp_tensor)
    pred_image_tensor = tf.reshape(pred_image, (1, width, height, 1));
    inv_warped_real_image = tf.reshape(batch_displacement_warp2d(pred_image_tensor, inverse_warp_tensor), pred_image_tensor.shape)

    sx = dx ;
    sy = dy
    for s in range(smooth_warp_iter):
        sx1 = smooth_image(sx)
        sy1 = smooth_image(sy)
        sx = sx1
        sy = sy1

    warp_tensor = tf.reshape(tf.stack([sx, sy]), (1, 2, width, height))
    warped_image_tensor = tf.reshape(batch_displacement_warp2d(pred_image_tensor, warp_tensor), real_image.shape)

    if image_weight > 0:
        dif = warped_image_tensor - real_image
        dif = dif ** 2
        image_loss = tf.sqrt(tf.reduce_mean(dif));
        dif = inv_warped_real_image - pred_image_tensor
        dif = dif ** 2
        image_loss += tf.sqrt(tf.reduce_mean(dif));
    else:
        image_loss = 0

    if sb_weight > 0:
        sb_loss = SBRegCF(real_image, warped_image_tensor) ;
        sb_loss += SBRegCF(inv_warped_real_image, pred_image_tensor) ;
    else:
        sb_loss = 0

    if (mag_weight > 0):
        mag_loss = compute_magnitude_loss(warp_tensor) ;
    else:
        mag_loss = 0

    if smooth_weight > 0:
        sx = smooth_image_niter(dx, smooth_iter)
        sy = smooth_image_niter(dy, smooth_iter)
        smooth_loss = tf.reduce_mean(tf.abs(tf.subtract(sx, dx)) + tf.abs(tf.subtract(sy, dy)));
    else:
        smooth_loss = 0

    if laplace_weight > 0:
        lx = laplace(dx);
        ly = laplace(dy);
        laplace_loss = tf.sqrt(tf.reduce_mean((lx + ly) ** 2))
    else:
        laplace_loss = 0

    if mi_weight > 0:
        Ixp, Iyp = sobel_filter2D(real_image) ;
        Ixw, Iyw = sobel_filter2D(warped_image_tensor) ;
        Iwmag = (Ixw ** 2 + Iyw ** 2) ;
        Ipmag = (Ixp**2 + Iyp**2) ;
        dif = abs((Ipmag - Iwmag)) 
        mi_loss = tf.sqrt(tf.reduce_mean(dif));
    else:
        mi_loss = 0

    jac = warp_jacobian(warp_tensor) 
    zoffset = 1e-9*tf.cast(jac<=0,tf.float32)

    negj = tf.cast(jac<=0,tf.float32) * jac
    negj *= -1 
    nzoffset = 1e-9*tf.cast(negj<=0,tf.float32)
    logj = tf.log((tf.cast(jac>0,tf.float32) * jac)+zoffset)
    neg_logj = tf.log((tf.cast(negj>0,tf.float32) * negj)+nzoffset)
    jac_loss = tf.reduce_mean(tf.abs(logj)) 
    NEG_AREA_K = 10.0
    NEG_AREA_K = 1.0
    MAX_NEG_RATIO = (40 / NEG_AREA_K)
    jac = tf.clip_by_value(jac, -MAX_NEG_RATIO, MAX_NEG_RATIO) ;
    jac_loss =  tf.reduce_mean((tf.log(1.0 + tf.exp(NEG_AREA_K * jac)) / NEG_AREA_K) - jac);

    return (jac_weight * jac_loss + mi_weight*mi_loss + image_weight*image_loss + mag_weight*mag_loss + laplace_weight*laplace_loss + smooth_weight*smooth_loss + sb_weight * sb_loss);


def make_kernel(a):
    """Transform a 2D array into a convolution kernel"""
    a = np.asarray(a)
    a = a.reshape(list(a.shape) + [1, 1])
    return tf.constant(a, dtype=1)


def simple_conv(x, k):
    """A simplified 2D convolution operation"""
    x = tf.expand_dims(tf.expand_dims(x, 0), -1)
    y = tf.nn.depthwise_conv2d(x, k, [1, 1, 1, 1], padding='SAME')
    return y[0, :, :, 0]


def laplace(x):
    """Compute the 2D laplacian of an array"""
    laplace_k = make_kernel([[0.5, 1.0, 0.5],
                             [1.0, -6., 1.0],
                             [0.5, 1.0, 0.5]])
    return simple_conv(x, laplace_k)


def smooth_image_niter(x, n):
    smooth_k = make_kernel([[1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0],
                            [1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0],
                            [1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0]])
    sx = x ;
    for iter in range(n):
        sx = simple_conv(sx, smooth_k)

    return sx

def smooth_image(x):
    """Compute the 3x3 average of an array"""
    smooth_k = make_kernel([[1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0],
                            [1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0],
                            [1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0]])
    return simple_conv(x, smooth_k)


def write_tensor_nii(t, fname):
    """write a tensor to a nifti volume"""
    if len(t.shape) > 3:
        if K.eval(t).shape[1] == 2:
            t = tf.transpose(t, [2, 3, 0, 1])
        elif K.eval(t).shape[3] != 3:  # don't squeeze a vector
            t = tf.squeeze(t)
    im = nib.Nifti1Image(K.eval(t), np.diag([1, 1, 1, 1]));
    nib.save(im, fname)


def loss_fn(y_true, y_pred):
    return y_true  # or y_pred


class LossHistory(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.warps = []

    def on_batch_end(self, batch, logs={}):
        self.losses.append(logs.get('loss'))


import keras.callbacks as cbks


class CustomMetrics(cbks.Callback):
    def on_epoch_end(self, epoch, logs=None):
        for k in logs:
            if k.endswith('warp_image_loss'):
                print(logs[k])

import numpy as np
def mutual_information(im1, im2, nbins):
     """ Mutual information for joint histogram
     """
     # Convert bins counts to probability values
     hgram = np.histogram2d(im1, im2, bins=nbins)
     pxy = hgram / float(np.sum(hgram))
     px = np.sum(pxy, axis=1) # marginal for x over y
     py = np.sum(pxy, axis=0) # marginal for y over x
     px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
     # Now we can do the calculation using the pxy, px_py 2D arrays
     nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
     return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

# def mutual_information(im1, im2, nbins):
#      """ Mutual information for joint histogram
#      """
#      # Convert bins counts to probability values
#      hgram = np.histogram2d(im1, im2, bins=nbins)
#      pxy = hgram / float(np.sum(hgram))
#      px = np.sum(pxy, axis=1) # marginal for x over y
#      py = np.sum(pxy, axis=0) # marginal for y over x
#      px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
#      # Now we can do the calculation using the pxy, px_py 2D arrays
#      nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
#      return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

def sobel_filter2D(im):
    """ Compute the sobel filter of the input image
    """
    sobel_x = make_kernel([[-1.0, 0, 1.0],
                           [-2.0, 0, 2.0],
                           [-1.0, 0, 1.0]])
    sobel_y = make_kernel([[-1.0, -2.0, -1.0],
                           [ 0.0,  0,    0.0],
                           [ 1.0,  2.0,  1.0]])
    im =  tf.cast(im, dtype=tf.float32)
    im_x = simple_conv(im, sobel_x)
    im_y = simple_conv(im, sobel_y)
    
    return im_x, im_y


def warp_inverse(warp_tensor, pad=15):
    """
    compute the inverse of the input displacement warp.
    """
    width = 512
    height = 512
    dx = warp_tensor[:,0,:,:]
    dy = warp_tensor[:,1,:,:]
    dxw =  tf.reshape(batch_displacement_warp2d(tf.reshape(dx, (1, width, height, 1)), warp_tensor), (width, height))
    dyw =  tf.reshape(batch_displacement_warp2d(tf.reshape(dy, (1, width, height, 1)), warp_tensor), (width, height))

#    x0 = tf.linspace(0.0, width-1.0, width)
#    y0 = tf.linspace(0.0, height-1.0, height)
#    p0x, p0y = tf.meshgrid(x0, y0)
#    p0xw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0x, (1, width, height, 1)), warp_tensor), (width, height))
#    p0yw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0y, (1, width, height, 1)), warp_tensor), (width, height))
#    dx = -2.*(p0xw-p0x) / (tf.cast(width, tf.float32))
#    dy = -2.*(p0yw-p0y) / (tf.cast(height,tf.float32))

    inverse_warp_tensor = tf.reshape(tf.stack([-1.*dxw, -1.*dyw]), (1, 2, width, height))
    return inverse_warp_tensor

def warp_jacobian(warp_tensor, pad=15):
    """
    compute the jacobian of the input warp.
    pad is the region on the boundary to 0
    """

    width = 512
    height = 512
    x0 = tf.linspace(0.0, width-1.0, width)
    y0 = tf.linspace(0.0, height-1.0, height)
    p0x, p0y = tf.meshgrid(x0, y0)
    
    p0xw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0x, (1, width, height, 1)), warp_tensor), (width, height))
    p0yw =  tf.reshape(batch_displacement_warp2d(tf.reshape(p0y, (1, width, height, 1)), warp_tensor), (width, height))
    
    p1xw = tf.manip.roll(p0xw, shift=[1,0], axis=[0,1]) ;
    p1yw = tf.manip.roll(p0yw, shift=[1,0], axis=[0,1]) ;
    
    p2xw = tf.manip.roll(p0xw, shift=[0,1], axis=[0,1]) ;
    p2yw = tf.manip.roll(p0yw, shift=[0,1], axis=[0,1]) ;
    
    de1x = p1xw - p0xw 
    de1y = p1yw - p0yw
    
    de2x = p2xw - p0xw 
    de2y = p2yw - p0yw

    jac = de2x*de1y - de2y*de1x
    pad0 = tf.cast(p0x > pad, tf.float32)
    pad1 = tf.cast(p0x < width-pad, tf.float32)
    pad2 = tf.cast(p0y > pad, tf.float32)
    pad3 = tf.cast(p0y < height-pad, tf.float32)
    pad = pad0*pad1*pad2*pad3
    return jac*pad
