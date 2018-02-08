from __future__ import print_function

from keras.models import Model
from keras.layers import Input, concatenate, \
    Conv2D, MaxPooling2D, Conv2DTranspose, Conv3D,MaxPooling3D, Activation,\
    Deconvolution3D,UpSampling3D, BatchNormalization, ZeroPadding3D, Cropping3D, MaxPool3D
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
    model = class_net(feature_shape=(256,256), dim=2, unet_num_filters=32, unet_depth=4, n_labels=2,
                      initial_learning_rate=0.001, loss='categorical_crossentropy')
    # model = atrous_net(input_shape=(16, 16, 16), num_filters=32)
    # model = unet_model_3d(input_shape = (16,16,16), num_filters=32, unet_depth=2,
    #                       downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=10,
    #                       loss='dice_coef_loss', initial_learning_rate=0.00001, deconvolution=False)

