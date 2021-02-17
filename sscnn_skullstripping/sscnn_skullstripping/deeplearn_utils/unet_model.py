from __future__ import print_function

from keras.models import Model
from keras.layers import Input, concatenate, \
    Conv2D, MaxPooling2D, Conv2DTranspose, Conv3D,MaxPooling3D, Activation,\
    Deconvolution3D,UpSampling3D, BatchNormalization, ZeroPadding3D, Cropping3D, MaxPool3D, AveragePooling3D
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint
from keras import backend as K
from keras.losses import  mean_squared_error, mean_absolute_percentage_error
from keras.models import Sequential
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten, UpSampling2D, Add
import keras
from keras.utils import  multi_gpu_model
import numpy as np
import os
import tensorflow as tf

from keras.utils.generic_utils import get_custom_objects

K.set_image_data_format('channels_last')


def dice_coef_loss2(y_true, y_pred):
    area_reg = 0.1

    y_true /= K.sum(y_true, axis=-1, keepdims=True)
    y_true = K.clip(y_true, K.epsilon(), 1)
    y_true_reshape = K.batch_flatten(y_true) #K.reshape(y_true, (num_samples, num_voxels, num_labels))

    y_pred /= K.sum(y_pred, axis=-1, keepdims=True)
    y_pred = K.clip(y_pred, K.epsilon(), 1)
    y_pred_reshape = K.batch_flatten(y_pred)#K.reshape(y_pred, (num_samples, num_voxels, num_labels))


    sum_over_axis = 1
    numerator = 2 * K.sum(y_true_reshape * y_pred_reshape, sum_over_axis)
    denominator = K.sum(K.square(y_true_reshape), sum_over_axis) + K.sum(K.square(y_pred_reshape), sum_over_axis)
    denominator = K.maximum(denominator, area_reg)

    dice_metric = numerator / denominator
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

    # input_shape_list = list(input_shape)
    # input_shape_list.append(1)
    # input_shape_append = tuple(input_shape_list)
    # print(input_shape)
    input_img = Input(shape=input_shape, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []
    # print('unet depth is ')
    # print(unet_depth)
    for i in range(unet_depth):

        prev = input_img if i == 0 else pools[i-1]
        # print(int(num_filters*(2**i)/downsize_filters_factor))
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
#     print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    # print(loss)
    # print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        # print(loss)
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
        # print('segmentation network')
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
    # print(input_shape_append)
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
    # print(input_shape_append)
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
    # print(input_shape_append)
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
    if dim == 2:
        ConvL =  Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2,2)
        UpSamplingL = UpSampling2D
        filter_shape = (5,5)
        out_filter_shape = (1,1)

    elif dim==3:
        ConvL = Conv3D
        MaxPoolingL= MaxPooling3D
        pool_size = (2,2,2)
        UpSamplingL = UpSampling3D
        filter_shape = (3,3,3)
        out_filter_shape = (1,1,1)

    # print('out filter shape is ' + str(out_filter_shape))
    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    # print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input' )
    convs = []
    pools = []
    inputs = []
    centered_inputs = []
    endpoints = []

    # print('unet depth is ')
    # print(unet_depth)
    for i in range(unet_depth):

        prev = input_img if i == 0 else pools[i-1]
        # print(int(num_filters*(2**i)/downsize_filters_factor))
        conv = ConvL(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_1_%d' % (i)))(prev)
        conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
        conv = ConvL(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_2_%d' % (i)))(conv)
        conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
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
        conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
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
    # print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    # print(loss)
    # print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        conv = ConvL(num_outputs, out_filter_shape, activation='linear',  name='final_conv_3d')(up)


        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=inputs, outputs=conv)
                parallel_model = multi_gpu_model(model, gpus=num_gpus)
                if loss == 'grad_loss':
                    # print(loss)
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
                else:
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)

                return model, parallel_model
        else:
            model = Model(inputs=inputs, outputs=conv)
            if loss == 'grad_loss':
                # print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        # print('segmentation network')
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






def unet_model_2d_noBN(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2), n_labels=0,
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
        # conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
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
    # print(convs)
    endpoints.append(convs[-1])


    up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
    print(loss)
    print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        conv = ConvL(num_outputs, out_filter_shape, activation='linear',  name='final_conv_3d')(up)


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






def conv_block(input_tensor, filters, stage, block, strides, dilation_rate):
    """A block that has a conv layer at shortcut.
    # Arguments
        input_tensor: input tensor
        kernel_size: default 3, the kernel size of middle conv layer at main path
        filters: list of integers, the filters of 3 conv layer at main path
        stage: integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names
        strides: Strides for the first conv layer in the block.
    # Returns
        Output tensor for the block.
    Note that from stage 3,
    the first conv layer at main path is with strides=(2, 2)
    And the shortcut should have strides=(2, 2) as well
    """

    dim = len(input_tensor.shape)
    #(num_samples, 32, 32,1)

    print(dim)
    if dim == 4:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        onexone_filter_shape = (1, 1)
        out_filter_shape = (1, 1)
        strides = strides

    #(num_samples, 32, 32, 32,1)
    elif dim == 5:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)
        strides = strides



    filters1, filters2, filters3 = filters
    if K.image_data_format() == 'channels_last':
        bn_axis = -1
    else:
        bn_axis = 1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = ConvL(filters1, onexone_filter_shape, padding='same', dilation_rate=dilation_rate,
               name=conv_name_base + '2a')(input_tensor)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = Activation('relu')(x)

    x = ConvL(filters2, filter_shape, padding='same',dilation_rate=dilation_rate,
               name=conv_name_base + '2b')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = Activation('relu')(x)

    x = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate, name=conv_name_base + '2c')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    shortcut = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate,
                      name=conv_name_base + '1')(input_tensor)
    shortcut = BatchNormalization(axis=bn_axis, name=bn_name_base + '1')(shortcut)

    x = Add()([shortcut, x])

    x = Activation('relu')(x)
    return x


def identity_block(input_tensor, filters, stage, block, dilation_rate):
    """The identity block is the block that has no conv layer at shortcut.
    # Arguments
        input_tensor: input tensor
        kernel_size: default 3, the kernel size of middle conv layer at main path
        filters: list of integers, the filters of 3 conv layer at main path
        stage: integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names
    # Returns
        Output tensor for the block.
    """

    dim = len(input_tensor.shape)
    if dim == 4:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        onexone_filter_shape = (1, 1)
        out_filter_shape = (1, 1)
        strides = (2,2)

    elif dim == 5:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)
        strides = (2,2,2)



    filters1, filters2, filters3 = filters
    if K.image_data_format() == 'channels_last':
        bn_axis = -1
    else:
        bn_axis = 1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = ConvL(filters1, onexone_filter_shape, padding='same',
              dilation_rate=dilation_rate, name=conv_name_base + '2a')(input_tensor)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = Activation('relu')(x)

    x = ConvL(filters2, filter_shape, padding='same',
              dilation_rate=dilation_rate, name=conv_name_base + '2b')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = Activation('relu')(x)

    x = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate,
              name=conv_name_base + '2c')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    x = Add()([x, input_tensor ])
    x = Activation('relu')(x)
    return x




def unet_model(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=0,
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

    print(input_shape)
    input_img = Input(shape=input_shape, name='input' )



    x = ConvL(32, (3, 3, 3), strides=(1,1,1), padding='same', name='conv1')(input_img)
    x = BatchNormalization(name='bn_conv1')(x)
    x = Activation('relu')(x)
    x = ConvL(32, (3, 3, 3), strides=(1,1,1), padding='same', name='conv2')(input_img)
    x = BatchNormalization(name='bn_conv2')(x)
    x = Activation('relu')(x)

    # x = MaxPooling3D((3, 3, 3), strides=(2, 2, 2))(x)
    #
    x = conv_block(x, [16, 16, 16], stage=2, block='a', strides=(1, 1, 1), dilation_rate=1)
    x = identity_block(x, [16, 16, 16], stage=2, block='b', dilation_rate=1)
    x = identity_block(x, [16, 16, 16], stage=2, block='c', dilation_rate=1)
    # #
    x = conv_block(x,  [32, 32, 32], stage=3, block='a', strides=(1, 1, 1), dilation_rate=2)
    x = identity_block(x, [32, 32, 32], stage=3, block='b', dilation_rate=2)
    # x = identity_block(x, [16, 32, 32], stagetage=3, block='c', dilation_rate=2)
    # x = identity_block(x, [16, 32, 32], stage=3, block='d', dilation_rate=2)
    # #
    x = conv_block(x, [16, 32, 32], stage=4, block='a', strides=(1, 1, 1), dilation_rate=4)
    x = identity_block(x,  [16, 32, 32], stage=4, block='b', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='c', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='d', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='e',  dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='f',  dilation_rate=4)

    x = conv_block(x,  [16, 32, 32], stage=5, block='a', strides=(1, 1, 1), dilation_rate=8)
    x = identity_block(x, [16, 32, 32], stage=5, block='b',  dilation_rate=8)
    # x = identity_block(x, [512, 512, 2048], stage=5, block='c',  dilation_rate=2)

    # decoding block


    print(loss)
    print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:
        print(loss)
        x = ConvL(1, (1, 1, 1), activation='relu', name='final_conv_3d')(x)

        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=input_img, outputs=x)
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
            model = Model(inputs=input_img, outputs=x)
            if loss == 'grad_loss':
                print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        print('segmentation network')
        if n_labels > 1:
            x = ConvL(n_labels, (1, 1, 1), activation='softmax', name='final_conv_3d')(x)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_img, outputs=x)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=input_img, outputs=x)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model
        else:
            x = ConvL(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(x)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_img, outputs=x)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=input_img, outputs=x)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model

    # x = AveragePooling3D((7, 7, 7), name='avg_pool')(x)
    # model = Model(input_img, x, name='resnet50')
    #
    # return model




#     print(input_shape)
#     input_img = Input(shape=input_shape, name='input' )
#     convs = []
#     pools = []
#     inputs = []
#     centered_inputs = []
#     endpoints = []
#     print('unet depth is ')
#     print(unet_depth)
#     for i in range(unet_depth):
#
#         prev = input_img if i == 0 else pools[i-1]
#         print(int(num_filters*(2**i)/downsize_filters_factor))
#         conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
#                       activation='relu', padding='same', kernel_initializer="he_normal",
#                       name=('conv3D_D_1_%d' % (i)))(prev)
#         conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
#         conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), filter_shape,
#                       activation='relu', padding='same', kernel_initializer="he_normal",
#                       name=('conv3D_D_2_%d' % (i)))(conv)
#         conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
#         if i < (unet_depth - 1):
#             pools.append(MaxPooling3D(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))
#
#         convs.append(conv)
#
#     for i in range(unet_depth - 1):
#         index = i + unet_depth - 1
#         level = unet_depth - (i + 2)
#         up = concatenate([UpSampling3D(size=pool_size,  name=('upsampling_U_%d' % (level+1)))(convs[index]),
#                           convs[level]], axis=-1,  name=('concat_%d' % (level)))
#         conv = Conv3D(num_filters * (2 ** level), (3, 3, 3), padding="same", activation="relu",
#                       kernel_initializer="he_normal",
#                       name=('conv3D_U_1_%d' % (level))
#                       )(up)
#         conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
#         conv = Conv3D(num_filters * (2 ** level), (3, 3, 3), padding="same", activation="relu",
#                       kernel_initializer="he_normal",
#                       name=('conv3D_U_2_%d' % (level)))(conv)
#         convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))
#
# #    conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
# #    conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu",
# #                  kernel_initializer="he_normal")(conv)
# #    conv = BatchNormalization()(conv)
# #    center_input = Cropping3D(cropping=(0, 0, 0))(input_img)
#
#     inputs.append(input_img)
# #    centered_inputs.append(center_input)
#     print(convs)
#     endpoints.append(convs[-1])
#
#
#     up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
#     print(loss)
#     print('is_seg_network' + str(is_seg_network))
#     if is_seg_network == False:
#         print(loss)
#         conv = Conv3D(1, (1,1,1), activation='relu',  name='final_conv_3d')(up)
#
#
#         if num_gpus > 1:
#             with tf.device('/cpu:0'):
#                 model = Model(inputs=inputs, outputs=conv)
#                 parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                 if loss == 'grad_loss':
#                     print(loss)
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
#                 else:
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
#
#                 return model, parallel_model
#         else:
#             model = Model(inputs=inputs, outputs=conv)
#             if loss == 'grad_loss':
#                 print(loss)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
#             else:
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
#             return model, model
#
#     else:
#         print('segmentation network')
#         if n_labels > 1:
#             conv = Conv3D(n_labels, (1, 1, 1), activation='softmax', name='final_conv_3d')(up)
#
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
#                                            loss=dice_coef_loss2, )
#
#                     model.compile(optimizer=Adam(lr=initial_learning_rate),
#                                   loss=dice_coef_loss2, )
#
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate),
#                               loss=dice_coef_loss2, )
#                 return model, model
#         else:
#             conv = Conv3D(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(up)
#
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
#                                            loss=dice_coef_loss2, )
#
#                     model.compile(optimizer=Adam(lr=initial_learning_rate),
#                                   loss=dice_coef_loss2, )
#
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate),
#                               loss=dice_coef_loss2, )
#                 return model, model





def build_oriented_unet_2d(input_layer, input_shape, num_filters, unet_depth, depth_per_level=1, layer_prefix='', downsize_filters_factor=1,
                         kernel_size=(3,3), pool_size=(2,2), n_labels=0, num_outputs=1, num_gpus=1):

    """

    :param input_img_ax: input layer
    :param input_shape: (256,256,1)
    :param num_filters: initial number of filters
    :param unet_depth: number of poolings
    :param downsize_filters_factor: generallly 1
    :param pool_size: (2,2)
    :param n_labels: number of labels to predict. 0 if regression
    :param num_outputs: dimensionality of outouts. 1 if single regression. 2 if vector valued output
    :return: callable final layer
    """
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
        filter_shape = kernel_size
        out_filter_shape = (1,1)

    elif dim==4:
        ConvL = Conv3D
        MaxPoolingL= MaxPooling3D
        pool_size = (2,2,2)
        UpSamplingL = UpSampling3D
        filter_shape = kernel_size
        out_filter_shape = (1,1,1)

    # print('out filter shape is ' + str(out_filter_shape))

    convs = []
    pools = []
    inputs = []
    endpoints = []

    # print('unet depth is ')
    # print(unet_depth)
    for i in range(unet_depth):

        prev = input_layer if i == 0 else pools[i - 1]
        # print(int(num_filters * (2 ** i) / downsize_filters_factor))
        conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
                     activation='relu', padding='same', kernel_initializer="he_normal",
                     name=('conv3D_D_1_%d' % (i) + str(0) + layer_prefix))(prev)
        conv = BatchNormalization(name=('bnorm_D_1_%d' % (i) + str(0) + layer_prefix))(conv)
        convs.append(conv)

        for iter in range(1,depth_per_level):

            conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
                         activation='relu', padding='same', kernel_initializer="he_normal",
                         name=('conv3D_D_1_%d' % (i) + str(iter) + layer_prefix))(convs[-1])
            conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)+ str(iter) + layer_prefix))(conv)
            convs.append(conv)

        if i < (unet_depth - 1):
            pools.append(MaxPoolingL(pool_size, name=('pool_D_%d' % (i)+ layer_prefix), data_format='channels_last')(conv))

            # convs.append(conv)
    # print(convs)
    # print(len(convs))
    for i in range(unet_depth - 1):
        # print(i)
        index = i*depth_per_level + unet_depth*depth_per_level - 1
        # print(index)
        level = unet_depth*depth_per_level - (i + 2)*depth_per_level + (depth_per_level-1)
        # print(level)
        # print(convs[index])
        # print(convs[level])
        actual_level = unet_depth - (i+2)
        up = concatenate([UpSamplingL(size=pool_size, name=('upsampling_U_%d' % (level + 1)+ layer_prefix))(convs[index]),
                          convs[level]], axis=-1, name=('concat_%d' % (level)+ layer_prefix))
        conv = ConvL(num_filters * (2 ** actual_level), filter_shape, padding="same", activation="relu",
                     kernel_initializer="he_normal",
                     name=('conv3D_U_1_%d' % (level) + str(0) + layer_prefix)
                     )(up)
        conv = BatchNormalization(name=('bnorm_U_1_%d' % (level) + str(0) + layer_prefix))(conv)
        convs.append(conv)

        for iter in range(1, depth_per_level):

            conv = ConvL(num_filters * (2 ** actual_level), filter_shape, padding="same", activation="relu",
                         kernel_initializer="he_normal",
                         name=('conv3D_U_2_%d' % (level) + str(iter) +  layer_prefix))(convs[-1])
            conv = BatchNormalization(name=('bnorm_U_2_%d' % (level)+ str(iter) + layer_prefix))(conv)
            convs.append(conv)

        #    conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
        #    conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu",
        #                  kernel_initializer="he_normal")(conv)
        #    conv = BatchNormalization()(conv)
        #    center_input = Cropping3D(cropping=(0, 0, 0))(input_img)

    inputs.append(input_layer)
    #    centered_inputs.append(center_input)
    # print(convs)
    endpoints.append(convs[-1])

    up = concatenate(inputs + endpoints, axis=-1, name='final_concat'+layer_prefix)

    # print('is_seg_network' + str(is_seg_network))
    if is_seg_network == False:

        conv = ConvL(num_outputs, out_filter_shape, activation='linear', name='final_conv_3d'+ layer_prefix)(up)

    else:
        # print('segmentation network')
        if n_labels > 1:
            conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d'+ layer_prefix)(up)
        else:
            conv = ConvL(1, out_filter_shape, activation='sigmoid', name='final_conv_3d'+ layer_prefix)(up)

    return conv


def build_compile_model(input_layer, input_shape, conv, n_labels,  loss, num_gpus, initial_learning_rate):
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


    if is_seg_network == False:
        # print(loss)
        # conv = ConvL(1, out_filter_shape, activation='relu',  name='final_conv_3d')(conv_penultimate)


        if num_gpus > 1:
            with tf.device('/cpu:0'):
                model = Model(inputs=input_layer, outputs=conv)
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
            model = Model(inputs=input_layer, outputs=conv)
            if loss == 'grad_loss':
                print(loss)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
            else:
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
            return model, model

    else:
        # print('segmentation network')
        if n_labels > 1:
            # conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d')(conv_penultimate)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_layer, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=input_layer, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model
        else:
            # conv = ConvL(1,out_filter_shape, activation='sigmoid', name='final_conv_3d')(conv_penultimate)

            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_layer, outputs=conv)

                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate),
                                           loss=dice_coef_loss2, )

                    model.compile(optimizer=Adam(lr=initial_learning_rate),
                                  loss=dice_coef_loss2, )

                    return model, parallel_model
            else:
                model = Model(inputs=input_layer, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate),
                              loss=dice_coef_loss2, )
                return model, model



def unet_2d_v1(input_shape, num_filters, unet_depth, depth_per_level=1, downsize_filters_factor=1, kernel_size=(3,3),pool_size=(2, 2), n_labels=0,
                  loss='mean_squared_error', initial_learning_rate=0.00001, deconvolution=False, num_gpus=1, num_outputs=1, add_modality_channel=False):
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    dim = len(input_shape)
    if add_modality_channel == True:
        #add 1 to the number of channels
        num_channels = input_shape[-1]
        linput_shape = list(input_shape)
        lshape = linput_shape[0:-1] + [linput_shape[-1]+2]
        input_shape = tuple(lshape)

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


    input_img_ax = Input(shape=input_shape, name='input_ax' )
    conv_ax = build_oriented_unet_2d(input_img_ax, input_shape=input_shape, num_filters=num_filters,
                                     unet_depth=unet_depth, depth_per_level=depth_per_level, layer_prefix='', downsize_filters_factor=downsize_filters_factor,
                                     kernel_size=kernel_size, pool_size=pool_size, n_labels=n_labels, num_outputs=num_outputs)
    (model_ax, parallel_model_ax) = build_compile_model(input_img_ax, input_shape,
                                                       conv_ax, n_labels,
                                                       loss, num_gpus, initial_learning_rate)
    return model_ax, parallel_model_ax




if __name__=="__main__":
    os.environ['LD_LIBRARY_PATH'] = '/usr/pubsw/packages/CUDA/lib64 '


    # model = unet_model(input_shape=(32,32,32,1), num_filters=32, unet_depth=4,
    #                    downsize_filters_factor=1,  n_labels=0, loss='mean_absolute_error',
    #                    initial_learning_rate=0.00001,
    #                    deconvolution=False,
    #                    use_patches=True, num_gpus=1)


    import numpy as np
    import keras
    import keras.utils
    import tensorflow as tf
    import keras.backend as K


    ytrue = np.random.randint(0, 3, (5, 4, 4, 4, 1))
    ytrue_soft = keras.utils.to_categorical(ytrue, 3)

    ypred = np.random.randint(0, 3, (5, 4, 4, 4, 1))
    ypred_soft = keras.utils.to_categorical(ypred, 3)

    tf_true = tf.constant(ytrue_soft, dtype=tf.float32)
    tf_pred = tf.constant(ypred_soft, dtype=tf.float32)




    def dice_coef_loss2(y_true, y_pred):
        area_reg = 0.1
        y_true /= K.sum(y_true, axis=-1, keepdims=True)
        y_true = K.clip(y_true, K.epsilon(), 1)

        # print(y_true.shape)


        batch_shape = y_true.shape.as_list()
        num_samples = batch_shape[0]

        num_labels = batch_shape[-1]

        num_voxels = np.prod(batch_shape[1:-1])
        y_true_reshape = K.reshape(y_true, (num_samples, num_voxels, num_labels))
        print(y_true_reshape.shape)


        y_pred /= K.sum(y_pred, axis=-1, keepdims=True)
        y_pred = K.clip(y_pred, K.epsilon(), 1)


        y_pred_reshape = K.reshape(y_pred, (num_samples, num_voxels, num_labels))
        print(y_pred_reshape.shape)


        sum_over_axis = 1
        numerator = 2 * K.sum(y_true_reshape * y_pred_reshape, sum_over_axis)
        print(numerator.shape)

        denominator = K.sum(K.square(y_true_reshape), sum_over_axis) + K.sum(K.square(y_pred_reshape), sum_over_axis)
        print(denominator.shape)
        # make sure we have no 0s on the bottom. K.epsilon()
        bottom = K.maximum(denominator, area_reg)
        print(bottom.shape)


        dice_metric = numerator / denominator
        print(dice_metric)
        dice_loss = 1 - dice_metric
        mean_dice_loss = K.mean(dice_loss)
        return mean_dice_loss


    dice = dice_coef_loss2(tf_true, tf_pred)
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        sess.run(dice)
        print(dice.eval())







    # model = atrous_net(input_shape=(16, 16, 16), num_filters=32)
    # model = unet_model_3d(input_shape = (16,16,16), num_filters=32, unet_depth=2,
    #                       downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=10,
    #                       loss='dice_coef_loss', initial_learning_rate=0.00001, deconvolution=False)

