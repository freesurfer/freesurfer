from __future__ import print_function

from keras.models import Model
from keras.layers import Input, concatenate, \
    Conv2D, MaxPooling2D, Conv2DTranspose, Conv3D,MaxPooling3D, Activation,\
    Deconvolution3D,UpSampling3D, BatchNormalization, ZeroPadding3D, Cropping3D
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint
from keras import backend as K

import numpy as np



K.set_image_data_format('channels_last')

def unet_model_3d(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=1,
                  initial_learning_rate=0.00001, deconvolution=False):
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

    for i in range(unet_depth):
        prev = input_img if i == 0 else pools[i-1]
        conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_1_%d' % (i)))(prev)
        conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
        conv = Conv3D(int(num_filters*(2**i)/downsize_filters_factor), (3, 3, 3),
                      activation='relu', padding='same', kernel_initializer="he_normal",
                      name=('conv3D_D_2_%d' % (i)))(conv)
        conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
        if i < (unet_depth - 1):
            pools.append(MaxPooling3D(pool_size, name=('pool_D_%d' % (i)))(conv))

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
#    endpoints.append(conv)

#    up = concatenate(centered_inputs + endpoints, axis=-1)
    conv = Conv3D(1, (1,1,1), activation='relu', name='final_conv_3d')(convs[-1])
    model = Model(inputs=inputs, outputs=conv)
    #
    #
    # conv1 = Conv3D(int(32/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(inputs)
    # conv1 = Conv3D(int(64/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv1)
    # pool1 = MaxPooling3D(pool_size=pool_size)(conv1)
    #
    # conv2 = Conv3D(int(64/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(pool1)
    # conv2 = Conv3D(int(128/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv2)
    # pool2 = MaxPooling3D(pool_size=pool_size)(conv2)
    #
    # conv3 = Conv3D(int(128/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(pool2)
    # conv3 = Conv3D(int(256/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv3)
    # pool3 = MaxPooling3D(pool_size=pool_size)(conv3)
    #
    # conv4 = Conv3D(int(256/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(pool3)
    # conv4 = Conv3D(int(512/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv4)
    #
    # up5 = get_upconv(pool_size=pool_size, deconvolution=deconvolution, depth=2,
    #                  nb_filters=int(512/downsize_filters_factor), image_shape=input_shape[-3:])(conv4)
    # up5 = concatenate([up5, conv3], axis=-1)
    # conv5 = Conv3D(int(256/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(up5)
    # conv5 = Conv3D(int(256/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv5)
    #
    # up6 = get_upconv(pool_size=pool_size, deconvolution=deconvolution, depth=1,
    #                  nb_filters=int(256/downsize_filters_factor), image_shape=input_shape[-3:])(conv5)
    # up6 = concatenate([up6, conv2], axis=-1)
    # conv6 = Conv3D(int(128/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(up6)
    # conv6 = Conv3D(int(128/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv6)
    #
    # up7 = get_upconv(pool_size=pool_size, deconvolution=deconvolution, depth=0,
    #                  nb_filters=int(128/downsize_filters_factor), image_shape=input_shape[-3:])(conv6)
    # up7 = concatenate([up7, conv1], axis=-1)
    # conv7 = Conv3D(int(64/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(up7)
    # conv7 = Conv3D(int(64/downsize_filters_factor), (3, 3, 3), activation='relu', padding='same')(conv7)
    #
    # conv8 = Conv3D(1, (1, 1, 1), activation="relu", name='final_conv_3d')(conv7)
    # model = Model(inputs, conv8)
    #conv8 = Conv3D(n_labels, (1, 1, 1))(conv7)
    #act = Activation('sigmoid')(conv8) # sigmoid for segmentation
    #model = Model(inputs=inputs, outputs=act)
    #model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss, metrics=[dice_coef])
    model.compile(optimizer=Adam(lr=initial_learning_rate), loss='mean_absolute_error')

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
    model = unet_model_3d((16,16,16), 32, 5, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=1,
                  initial_learning_rate=0.00001, deconvolution=False)
