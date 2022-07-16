import os
import numpy as np
import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import *
from tensorflow.keras.losses import mse
#from tensorflow.keras.utils import multi_gpu_model
from tensorflow.python.keras.utils.multi_gpu_utils import multi_gpu_model
from ._utility import dice_coef_loss2, grad_loss
import pdb as gdb

K.set_image_data_format('channels_last')


def sampling(args):
    """Reparameterization trick by sampling fr an isotropic unit Gaussian.

    # Arguments:
        args (tensor): mean and log of variance of Q(z|X)

    # Returns:
        z (tensor): sampled latent vector
    """

    z_mean, z_log_var = args
    batch = K.shape(z_mean)[0]
    dim = K.int_shape(z_mean)[1]
    # by default, random_normal has mean=0 and std=1.0
    epsilon = K.random_normal(shape=(batch, dim))
    return z_mean + K.exp(0.5 * z_log_var) * epsilon


def vae_encoder(input_shape, num_filters, num_poolings, latent_dim):

    dim = len(input_shape)
    num_channels = input_shape[-1]

    if dim == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        out_filter_shape = (1, 1)
        onexone_filter_shape = (1, 1)

    elif dim == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        out_filter_shape = (1, 1, 1)
        onexone_filter_shape = (1, 1, 1)

    elif dim == 2:
        ConvL = Conv1D
        MaxPoolingL = MaxPooling1D
        pool_size = (2)
        UpSamplingL = UpSampling1D
        filter_shape = (3)
        out_filter_shape = (1)
        onexone_filter_shape = (1)



    input_layer = Input(shape=input_shape, name='input')

    convs = []
    pools = []
    inputs = []
    endpoints = []

    print('num poolings are  ')
    print(num_poolings)

    for i in range(num_poolings):
        if i == 0:
            if num_channels > 100:
                conv1x1 = ConvL(int(num_filters * (2 ** 0)), onexone_filter_shape, padding='same', dilation_rate=1,
                                name='2a' + str(0))(
                    input_layer)
                conv1x1 = BatchNormalization()(conv1x1)
                conv1x1 = Activation('relu')(conv1x1)
                prev = conv1x1
            else:
                prev = input_layer
        else:
            prev = pools[i - 1]

        conv = ConvL(int(num_filters * (2 ** i) / 1), filter_shape,
                     activation='relu', padding='same', kernel_initializer="he_normal",
                     name=('conv3D_D_1_%d' % (i) ))(prev)
        # conv = BatchNormalization(name=('bnorm_D_1_%d' % (i) ))(conv)
        conv = ConvL(int(num_filters * (2 ** i) / 1), filter_shape,
                     activation='relu', padding='same', kernel_initializer="he_normal",
                     name=('conv3D_D_2_%d' % (i) ))(conv)
        # conv = BatchNormalization(name=('bnorm_D_2_%d' % (i) ))(conv)
        print(i)
        # conv = MaxPoolingL(pool_size, name=('pool_D_%d' % (i) ), data_format='channels_last')(conv)
        conv = MaxPoolingL(pool_size, name=('pool_D_%d' % (i) ), data_format='channels_last')(conv)
        pools.append(conv)
        convs.append(conv)

    # conv = pools[-1](conv)
    # add the
    shape = K.int_shape(conv)

    # generate latent vector Q(z|X)
    conv = Flatten()(conv)
    # conv = Dense(2*latent_dim, activation='relu', kernel_initializer="he_normal")(conv)
    z_mean = Dense(latent_dim, name='z_mean')(conv)
    z_log_var = Dense(latent_dim, name='z_log_var')(conv)

    # use reparameterization trick to push the sampling out as input
    # note that "output_shape" isn't necessary with the TensorFlow backend
    z = Lambda(sampling, output_shape=(latent_dim,), name='z')([z_mean, z_log_var])

    # instantiate encoder model
    encoder = Model(input_layer, [z_mean, z_log_var, z], name='encoder')
    pre_sample_shape = shape[1:]
    return encoder, pre_sample_shape


def vae_decoder(latent_dim, pre_sample_shape, num_filters, num_poolings, output_shape):
    dim = len(output_shape)
    num_channels = output_shape[-1]

    if dim == 3:
        ConvL = Conv2D
        ConvLTranspose = Conv2DTranspose
        MaxPoolingL = MaxPooling2D

        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        out_filter_shape = (1, 1)
        onexone_filter_shape = (1, 1)

    elif dim == 4:
        ConvL = Conv3D
        ConvLTranspose = Conv3DTranspose
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        out_filter_shape = (1, 1, 1)
        onexone_filter_shape = (1, 1, 1)

    elif dim == 2:
        ConvL = Conv1D
      # need something for conv1dtranspose...(maybe upsampling)
        MaxPoolingL = MaxPooling1D
        pool_size = (2)
        UpSamplingL = UpSampling1D
        filter_shape = (3)
        out_filter_shape = (1)
        onexone_filter_shape = (1)


    latent_inputs = Input(shape=(latent_dim,), name='z_sampling')
    flat_shape = np.prod(pre_sample_shape)

    conv = Dense(flat_shape, activation='relu')(latent_inputs)
    conv = Reshape((pre_sample_shape[0], pre_sample_shape[1], pre_sample_shape[2]))(conv)

    for i in range(num_poolings):
        level = num_poolings - i - 1
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                     kernel_initializer="he_normal",
                     name=('conv3D_U_1_%d' % (level)))(conv)
        # conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                     kernel_initializer="he_normal",
                     name=('conv3D_U_2_%d' %  (level)))(conv)
        # conv = BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv)
        conv = ConvLTranspose(filters=num_filters * (2 ** level),
                            kernel_size=filter_shape,
                            activation='relu',
                            strides=2,
                            padding='same')(conv)

    outputs = ConvLTranspose(filters=1,
                              kernel_size=filter_shape,
                              activation='relu',
                              padding='same',
                              name='decoder_output')(conv)
    decoder = Model(latent_inputs, outputs, name='decoder')
    return decoder



def vae(input_shape, num_filters, num_poolings, latent_dim, initial_learning_rate=0.001):

    encoder, pre_sample_shape = vae_encoder(input_shape=input_shape,  num_filters=num_filters,
                          num_poolings=num_poolings, latent_dim=latent_dim)

    decoder = vae_decoder(latent_dim=latent_dim, pre_sample_shape=pre_sample_shape, num_filters=num_filters,
                          num_poolings=num_poolings, output_shape=input_shape)

    outputs = decoder(encoder(encoder.inputs)[2])
    vae = Model(encoder.inputs, outputs, name='vae')

    reconstruction_loss = mse(K.flatten(encoder.inputs), K.flatten(outputs))
    reconstruction_loss *= np.prod(input_shape)
    kl_loss = 1 + encoder.get_layer('z_log_var').output - K.square(encoder.get_layer('z_mean').output) - K.exp(encoder.get_layer('z_log_var').output)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer=Adam(lr=initial_learning_rate), loss=None)
    vae.summary()

    return vae, encoder, decoder






def unet(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=None, n_labels=0,
               loss='mean_squared_error', initial_learning_rate=0.00001, deconvolution=False, num_gpus=1,
               num_outputs=1,GPnet=None,pooling='max', collapse_dim=None, batch_norm = True,
         output_shape = None):
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    assert unet_depth <= np.log2(np.array(input_shape[0:-1]).min()), 'unet too deep (%d) for min input patch shape %d' % (unet_depth, np.array(input_shape[0:-1]).min())
        
    dim = len(input_shape)
    if dim == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        if pool_size == None:
            pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        out_filter_shape = (1, 1)
    elif dim == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        if pool_size == None:
            pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        out_filter_shape = (1, 1, 1)

    input_img_ax = Input(shape=input_shape, name='input_ax')
    conv_ax = build_oriented_unet(input_img_ax, input_shape=input_shape, num_filters=num_filters,
                                     unet_depth=unet_depth, layer_prefix='',
                                     downsize_filters_factor=downsize_filters_factor,batch_norm=batch_norm,
                                  pool_size=pool_size, n_labels=n_labels, num_outputs=num_outputs,
                                  GPnet=GPnet, pooling=pooling,collapse_dim=collapse_dim,
                                  output_shape=output_shape)
    if (output_shape != None):
        vec_mse_loss = lambda yt, yp: K.sqrt(K.mean(K.sum(K.square(yt - yp), -1)))
        loss = [vec_mse_loss]

    (model_ax, parallel_model_ax) = build_compile_model(input_img_ax, input_shape, conv_ax, n_labels,
                                                        loss, num_gpus, initial_learning_rate)
    return model_ax, parallel_model_ax


def build_oriented_unet(input_layer, input_shape, num_filters, unet_depth, layer_prefix='',
                        downsize_filters_factor=1,
                        pool_size=None, n_labels=0, num_outputs=1, num_gpus=1, GPnet=None,pooling='max',
                        collapse_dim=None, batch_norm=True,output_shape=None):
    """
    Args:
        input_img_ax: input layer
        input_shape: (256,256,1)
        num_filters: initial number of filters
        unet_depth: number of poolings
        downsize_filters_factor: generally 1
        pool_size: (2,2)
        n_labels: number of labels to predict. 0 if regression
        num_outputs: dimensionality of outputs. 1 if single regression. 2 if vector valued output

    Returns:
        callable final layer

    """
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    dim = len(input_shape)
    num_channels = input_shape[-1]

    if dim == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        AvgPoolingL = AveragePooling2D
        if pool_size == None:
            pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (3, 3)
        onexone_filter_shape = (1, 1)
        out_filter_shape = (1, 1)
    elif dim == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        AvgPoolingL = AveragePooling3D
        if pool_size == None:
            pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)

    if pooling == 'max':
        PoolingL = MaxPoolingL
    else:
        PoolingL = AvgPoolingL
    print('out filter shape is ' + str(out_filter_shape))

    convs = []
    pools = []
    inputs = []
    endpoints = []

    print('unet depth is ')
    print(unet_depth)

    for i in range(unet_depth):
        if i == 0:
            if num_channels > 100000:
                conv1x1 = ConvL(int(num_filters * (2 ** 5)), onexone_filter_shape, padding='same', dilation_rate=1,
                                name='2a' + str(0))(
                    input_layer)
                if batch_norm == True:
                    conv1x1 = BatchNormalization()(conv1x1)
                conv1x1 = Activation('relu')(conv1x1)
                prev = conv1x1
            else:
                prev = input_layer
        else:
            prev = pools[i - 1]
        print(int(num_filters * (2 ** i) / downsize_filters_factor))

        conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
                     activation='relu', padding='same', kernel_initializer="he_normal",
                     name=('conv3D_D_1_%d' % (i) + layer_prefix))(prev)
        if batch_norm == True:
            conv = BatchNormalization(name=('bnorm_D_1_%d' % (i) + layer_prefix))(conv)
        conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
                     activation='relu', padding='same', kernel_initializer="he_normal",
                     name=('conv3D_D_2_%d' % (i) + layer_prefix))(conv)
        if batch_norm == True:
            conv = BatchNormalization(name=('bnorm_D_2_%d' % (i) + layer_prefix))(conv)
        if i < (unet_depth - 1):
            pools.append(
                PoolingL(pool_size, name=('pool_D_%d' % (i) + layer_prefix), data_format='channels_last')(conv))
        convs.append(conv)

    for i in range(unet_depth - 1):
        index = i + unet_depth - 1
        level = unet_depth - (i + 2)
        up = concatenate(
            [UpSamplingL(size=pool_size, name=('upsampling_U_%d' % (level + 1) + layer_prefix))(convs[index]),
             convs[level]], axis=-1, name=('concat_%d' % (level) + layer_prefix))
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                     kernel_initializer="he_normal",
                     name=('conv3D_U_1_%d' % (level) + layer_prefix))(up)
        if batch_norm == True:
            conv = BatchNormalization(name=('bnorm_U_1_%d' % (level) + layer_prefix))(conv)
        conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
                     kernel_initializer="he_normal",
                     name=('conv3D_U_2_%d' % (level) + layer_prefix))(conv)
        if batch_norm == True:
            convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level) + layer_prefix))(conv))
        else:
            convs.append(conv)

        # conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
        # conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu", kernel_initializer="he_normal")(conv)
        # conv = BatchNormalization()(conv)
        # center_input = Cropping3D(cropping=(0, 0, 0))(input_img)

    inputs.append(input_layer)
    # centered_inputs.append(center_input)
    print(convs)
    endpoints.append(convs[-1])
    upL = concatenate(inputs + endpoints, axis=-1, name='final_concat' + layer_prefix)

    if GPnet == 'max':
        upL = MaxPoolingL(conv.shape[1:dim], name=('pool_GP' + layer_prefix), data_format='channels_last')(upL)
    elif GPnet == 'avg':
        upL = AvgPoolingL(conv.shape[1:dim], name=('pool_GP' + layer_prefix), data_format='channels_last')(upL)

    if collapse_dim != None:
        upL = AvgPoolingL((1,1,conv.shape[dim-1]), name=('pool_GP' + layer_prefix), data_format='channels_last')(upL)
        ConvL = Conv2D
        upL_shape = K.int_shape(upL)
        upL = Reshape((upL_shape[1], upL_shape[2], upL_shape[4]))(upL)
        out_filter_shape = (1,1)
#        upL = ConvL(1, (1,1), padding="same", activation="relu", kernel_initializer="he_normal",name='conv3D_to_2D')(upL)
#        upL = Reshape(tuple([upL.shape[0:collapse_dim+1]])+tuple([upL.shape[dim]]))(upL)

        
    print('is_seg_network' + str(is_seg_network))
    use_bias = True
    if is_seg_network == False:
        conv = ConvL(1, out_filter_shape, activation='linear', name='final_conv_3d' + layer_prefix)(upL)
        if (output_shape != None):
            conv = AvgPoolingL((4,4,4), name=('pool_before_reshape' + layer_prefix), data_format='channels_last')(conv)
            conv = Flatten()(conv)
            ndown = 3
            oshape = tuple((int(output_shape[0]/2**ndown), int(output_shape[1]/2**ndown),output_shape[2]))
            conv = Dense(np.array(oshape).prod())(conv)
            conv = Reshape(oshape)(conv)
            upsample_shape = tuple((int(2**ndown),int(2**ndown),1))
            upsample_shape = (2,2,1)
            conv = Conv2D(num_filters, (3,3),
                     activation='relu', padding='same', kernel_initializer="he_normal",
                         name=('conv2D_reshape_0'))(conv)
            for d in range(ndown):
                conv = UpSampling2D(size=(2,2), name=('upsampling_after_reshape_%d' % d))(conv)

                conv = Conv2D(num_filters, (3,3),
                             activation='relu', padding='same', kernel_initializer="he_normal",
                             name=('conv2D_reshape_%d' % (d+1)))(conv)

            conv = Conv2D(output_shape[2], (1,1), activation='linear', name='final_conv_2d')(conv)

    else:
        print('segmentation network')
        if n_labels > 1:
            if (num_outputs > 1):
                ishape = K.int_shape(upL)+tuple([1])
                reL = ConvL(num_outputs*n_labels, out_filter_shape, activation='softmax', name='final_conv_3d' + layer_prefix)(upL)
                conv = Reshape((ishape[1], ishape[2], ishape[3], num_outputs,n_labels))(reL)
            else:
                conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d' + layer_prefix)(upL)
        else:
            if (GPnet != None):
                use_bias = False

            conv = ConvL(1, out_filter_shape, activation='sigmoid', name='final_conv_3d' + layer_prefix,use_bias=use_bias)(upL)

    return conv


def build_compile_model(input_layer, input_shape, conv, n_labels, loss, num_gpus, initial_learning_rate):
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

#    dim = len(input_shape)
#    if dim == 3:
#        ConvL = Conv2D
#        filter_shape = (5, 5)
#        out_filter_shape = (1, 1)
#
#    elif dim == 4:
#        ConvL = Conv3D
#        filter_shape = (3, 3, 3)
#        out_filter_shape = (1, 1, 1)

    if is_seg_network == False:
        print(loss)
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
        print('segmentation network')
        if n_labels > 1:
            # conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d')(conv_penultimate)
            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_layer, outputs=conv)
                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    return model, parallel_model
            else:
                model = Model(inputs=input_layer, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                return model, model
        else:
            # conv = ConvL(1,out_filter_shape, activation='sigmoid', name='final_conv_3d')(conv_penultimate)
            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_layer, outputs=conv)
                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    return model, parallel_model
            else:
                model = Model(inputs=input_layer, outputs=conv)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                return model, model


def unet_encoder_dense(feature_shape, output_shape, unet_num_filters, depth, depth_per_level,
                       n_labels, initial_learning_rate, loss='binary_crossentropy'):
    if len(feature_shape) == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        filter_shape = (3, 3)
        pool_shape = (2, 2)
        conv1x1_shape = (1, 1)
    elif len(feature_shape) == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        filter_shape = (3, 3, 3)
        pool_shape = (2, 2, 2)
        conv1x1_shape = (1, 1, 1)

    if n_labels != output_shape[-1]:
        print('Number of labels do not match the output shape last dimension')
        return None

    min_feature_dim = np.min(feature_shape[0:-1])
    if 2 ** depth > min_feature_dim:
        print('Reduce depth of the network to fit ' + str(depth) + ' poolings')
        return None

    num_output_dense_filters = np.prod(output_shape)

    input_shape_list = list(feature_shape)
    # input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    model = Sequential()
    for iter_layer in range(depth):
        if iter_layer == 0:
            model.add(ConvL(unet_num_filters * (2 ** iter_layer), filter_shape, padding='same', activation='relu',
                            input_shape=input_shape_append, kernel_initializer="he_normal"))
            for iter_depth_per_layer in range(depth_per_level - 1):
                model.add(ConvL(unet_num_filters * (2 ** iter_layer), filter_shape, padding='same',
                                activation='relu', kernel_initializer="he_normal"))
            model.add(MaxPoolingL(pool_size=pool_shape))
            model.add(Dropout(0.25))
        else:
            for iter_depth_per_layer in range(depth_per_level):
                model.add(ConvL(unet_num_filters * (2 ** iter_layer), filter_shape, padding='same',
                                activation='relu', kernel_initializer="he_normal"))
            model.add(MaxPoolingL(pool_size=pool_shape))
            model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(512, activation='relu', kernel_initializer="he_normal"))
    model.add(Dropout(0.5))
    model.add(Dense(256, activation='relu', kernel_initializer="he_normal"))
    model.add(Dropout(0.5))

    if n_labels == 1:
        model.add(Dense(num_output_dense_filters, activation='relu', kernel_initializer="he_normal"))
        model.add(Reshape(output_shape))
    elif n_labels > 1:
        model.add(Dense(num_output_dense_filters, activation='relu', kernel_initializer="he_normal"))
        model.add(Reshape(output_shape))
        model.add(ConvL(n_labels, conv1x1_shape, activation='softmax', kernel_initializer="he_normal"))

    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
    return model




def classnet(feature_shape, unet_num_filters, depth, depth_per_level, n_labels, initial_learning_rate,
                 loss='binary_crossentropy', batch_norm=True, dropout_frac=None,pooling='max',num_outputs=1,final_activation='softmax', conv_size=3, dense_size=512, feature_scale_with_depth=2):
    dim = len(feature_shape)
    num_channels = feature_shape[-1]
    if dim == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        AvgPoolingL = AveragePooling2D
        pool_shape = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (conv_size, conv_size)
        out_filter_shape = (1, 1)
    elif dim == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        AvgPoolingL = AveragePooling3D
        pool_shape = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (conv_size,conv_size,conv_size)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)
    elif dim == 2:
        out_filter_shape = (1)
        ConvL = Conv1D
        MaxPoolingL = MaxPooling1D
        onexone_filter_shape = (1, 1)
        AvgPoolingL = AveragePooling1D
        pool_shape = (2)
        UpSamplingL = UpSampling1D
        filter_shape = (conv_size)
        out_filter_shape = (1)

    if pooling == 'max':
        PoolingL = MaxPoolingL
    else:
        PoolingL = AvgPoolingL
    min_feature_dim = np.min(feature_shape[0:-1])
    if 2 ** depth > min_feature_dim:
        print('Reduce depth of the network to fit ' + str(depth) + ' poolings')
        return None

    input_shape_list = list(feature_shape)
    # input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    model = Sequential()
    use_bias = True

    # if very high number of channels, use a 1x1 filter to reduce the dimensionality
    if num_channels > 100 and 0:
        model.add(ConvL(int(unet_num_filters * (2 ** 0)), onexone_filter_shape, padding='same', dilation_rate=1,
                        name='2a' + str(0)))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

    unet_num_filters_at_depth = unet_num_filters
    for iter_layer in range(depth):
        if iter_layer == 0:

            model.add(ConvL(unet_num_filters_at_depth, filter_shape, padding='same', activation='relu',
                            input_shape=input_shape_append, kernel_initializer="he_normal", use_bias=use_bias))
            if batch_norm == True:
                model.add(BatchNormalization())

            for iter_depth_per_layer in range(depth_per_level - 1):
                model.add(ConvL(unet_num_filters_at_depth, filter_shape, padding='same', activation='relu',
                                kernel_initializer="he_normal", use_bias=use_bias))
                if batch_norm == True:
                    model.add(BatchNormalization())

            model.add(PoolingL(pool_size=pool_shape))
            if (dropout_frac != None):
                model.add(Dropout(dropout_frac))
        else:
            for iter_depth_per_layer in range(depth_per_level):
                model.add(ConvL(unet_num_filters_at_depth, filter_shape, padding='same', activation='relu',
                                kernel_initializer="he_normal", use_bias=use_bias))
                if batch_norm == True:
                    model.add(BatchNormalization())

            model.add(PoolingL(pool_size=pool_shape))
            if (dropout_frac != None):
                model.add(Dropout(dropout_frac))
        unet_num_filters_at_depth = int(unet_num_filters_at_depth*feature_scale_with_depth)

    model.add(Flatten())
    model.add(Dense(dense_size, activation='relu', kernel_initializer="he_normal", use_bias=use_bias))
    model.add(Dropout(0.2))
    if n_labels > 0:
        model.add(Dense(n_labels, activation=final_activation, use_bias=use_bias, kernel_initializer="he_normal"))
        model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss, metrics=['accuracy'])

    elif n_labels == 0:
        model.add(Dense(num_outputs, activation=final_activation, kernel_initializer="he_normal", use_bias=use_bias))
        model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)

    return model


#

def atrous_net(input_shape, num_filters, initial_learning_rate=0.00001, loss='mean_absolute_error'):
    input_shape_list = list(input_shape)
    input_shape_list.append(1)
    input_shape_append = tuple(input_shape_list)
    print(input_shape_append)
    input_img = Input(shape=input_shape_append, name='input')
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


def resnet_model(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=0,
                 loss='mean_absolute_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True,
                 num_gpus=1):
    if n_labels > 0:
        is_seg_network = True
    else:
        is_seg_network = False

    dim = len(input_shape)
    if dim == 3:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        out_filter_shape = (1, 1)
    elif dim == 4:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        out_filter_shape = (1, 1, 1)

    print(input_shape)
    input_img = Input(shape=input_shape, name='input')

    x = ConvL(16, (3, 3, 3), padding='same', name='conv1')(input_img)
    x = BatchNormalization(name='bn_conv1')(x)
    x = Activation('relu')(x)
    x = ConvL(16, (3, 3, 3), padding='same', name='conv2')(input_img)
    x = BatchNormalization(name='bn_conv2')(x)
    x = Activation('relu')(x)

    # x = niftynet_block(x, [16, 16], stage=2, block='b', dilation_rate=1)
    # x = niftynet_block(x, [16, 16], stage=2, block='c', dilation_rate=1)

    # x = niftynet_block(x, [32, 32], stage=3, block='b', dilation_rate=2)
    # x = niftynet_block(x, [32, 32], stage=3, block='c', dilation_rate=2)

    # x = niftynet_block(x, [64, 64], stage=4, block='b', dilation_rate=4)
    # x = niftynet_block(x, [64, 64], stage=4, block='c', dilation_rate=4)

    # x = MaxPooling3D((3, 3, 3), strides=(2, 2, 2))(x)

    x = conv_block(x, [16, 16, 16], stage=2, block='a', strides=(1, 1, 1), dilation_rate=1)
    x = conv_block(x, [16, 16, 16], stage=2, block='b', strides=(1, 1, 1), dilation_rate=1)
    x = conv_block(x, [16, 16, 16], stage=2, block='c', strides=(1, 1, 1), dilation_rate=1)

    # x = identity_block(x, [16, 16, 16], stage=2, block='b', dilation_rate=1)
    # x = identity_block(x, [16, 16, 16], stage=2, block='c', dilation_rate=1)

    x = conv_block(x, [32, 32, 32], stage=3, block='a', strides=(1, 1, 1), dilation_rate=2)
    x = conv_block(x, [32, 32, 32], stage=3, block='b', strides=(1, 1, 1), dilation_rate=2)
    x = conv_block(x, [32, 32, 32], stage=3, block='c', strides=(1, 1, 1), dilation_rate=2)

    x = conv_block(x, [64, 64, 64], stage=4, block='a', strides=(1, 1, 1), dilation_rate=4)
    x = conv_block(x, [64, 64, 64], stage=4, block='b', strides=(1, 1, 1), dilation_rate=4)
    x = conv_block(x, [64, 64, 64], stage=4, block='c', strides=(1, 1, 1), dilation_rate=4)

    # x = identity_block(x, [32, 32, 32], stage=3, block='b', dilation_rate=2)
    # x = identity_block(x, [32, 32, 32], stage=3, block='c', dilation_rate=2)
    # x = identity_block(x, [16, 32, 32], stage=3, block='d', dilation_rate=2)

    # x = conv_block(x, [16, 32, 32], stage=4, block='a', strides=(1, 1, 1), dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='b', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='c', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='d', dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='e',  dilation_rate=4)
    # x = identity_block(x,  [16, 32, 32], stage=4, block='f',  dilation_rate=4)

    # x = conv_block(x,  [16, 32, 32], stage=5, block='a', strides=(1, 1, 1), dilation_rate=8)
    # x = identity_block(x, [16, 32, 32], stage=5, block='b',  dilation_rate=8)
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
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    return model, parallel_model
            else:
                model = Model(inputs=input_img, outputs=x)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                return model, model
        else:
            x = ConvL(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(x)
            if num_gpus > 1:
                with tf.device('/cpu:0'):
                    model = Model(inputs=input_img, outputs=x)
                    parallel_model = multi_gpu_model(model, gpus=num_gpus)
                    parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                    return model, parallel_model
            else:
                model = Model(inputs=input_img, outputs=x)
                model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
                return model, model


def identity_block(input_tensor, filters, stage, block, dilation_rate):
    """The identity block is the block that has no conv layer at shortcut.

    Args:
        input_tensor: input tensor
        kernel_size: default 3, the kernel size of middle conv layer at main path
        filters: list of integers, the filters of 3 conv layer at main path
        stage: integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names

    Returns:
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
        strides = (2, 2)
    elif dim == 5:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)
        strides = (2, 2, 2)

    filters1, filters2, filters3 = filters
    if K.image_data_format() == 'channels_last':
        bn_axis = -1
    else:
        bn_axis = 1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = ConvL(filters1, onexone_filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2a')(
        input_tensor)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = Activation('relu')(x)

    x = ConvL(filters2, filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2b')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = Activation('relu')(x)

    x = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate, name=conv_name_base + '2c')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    x = Add()([x, input_tensor])
    x = Activation('relu')(x)
    return x


def conv_block(input_tensor, filters, stage, block, strides, dilation_rate):
    """A block that has a conv layer at shortcut.

    Note:
        From stage 3, the first conv layer at main path is with strides=(2, 2).
        The shortcut should have strides=(2, 2) as well

    Arguments:
        input_tensor: Input tensor
        filters: List of integers, the filters of 3 conv layer at main path
        stage: Integer, current stage label, used for generating layer names
        block: 'a','b'..., current block label, used for generating layer names
        strides: Strides for the first conv layer in the block

    Returns:
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
        strides = strides
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

    x = ConvL(filters1, onexone_filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2a')(
        input_tensor)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = Activation('relu')(x)

    x = ConvL(filters2, filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2b')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = Activation('relu')(x)

    x = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate, name=conv_name_base + '2c')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2c')(x)

    shortcut = ConvL(filters3, onexone_filter_shape, dilation_rate=dilation_rate, name=conv_name_base + '1')(
        input_tensor)
    shortcut = BatchNormalization(axis=bn_axis, name=bn_name_base + '1')(shortcut)

    x = Add()([shortcut, x])
    x = Activation('relu')(x)
    return x


def niftynet_block(input_tensor, filters, stage, block, dilation_rate):
    dim = len(input_tensor.shape)
    if dim == 4:
        ConvL = Conv2D
        MaxPoolingL = MaxPooling2D
        pool_size = (2, 2)
        UpSamplingL = UpSampling2D
        filter_shape = (5, 5)
        onexone_filter_shape = (1, 1)
        out_filter_shape = (1, 1)
        strides = (2, 2)
    elif dim == 5:
        ConvL = Conv3D
        MaxPoolingL = MaxPooling3D
        pool_size = (2, 2, 2)
        UpSamplingL = UpSampling3D
        filter_shape = (3, 3, 3)
        onexone_filter_shape = (1, 1, 1)
        out_filter_shape = (1, 1, 1)
        strides = (2, 2, 2)

    filters1, filters2 = filters
    if K.image_data_format() == 'channels_last':
        bn_axis = -1
    else:
        bn_axis = 1
    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    x = ConvL(filters1, filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2a')(
        input_tensor)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2a')(x)
    x = Activation('relu')(x)

    x = ConvL(filters2, filter_shape, padding='same', dilation_rate=dilation_rate, name=conv_name_base + '2b')(x)
    x = BatchNormalization(axis=bn_axis, name=bn_name_base + '2b')(x)
    x = Activation('relu')(x)

    x = Add()([x, input_tensor])
    x = Activation('relu')(x)
    return x


def compute_level_output_shape(filters, depth, pool_size, image_shape):
    """
    Each level has a particular output shape based on the number of filters used in that level and the depth or number
    of max pooling operations that have been done on the data at that point.

    Args:
        image_shape: shape of the 3d image.
        pool_size: the pool_size parameter used in the max pooling operation.
        filters: Number of filters used by the last node in a given level.
        depth: The number of levels down in the U-shaped model a given node is.

    Returns:
        A 5D vector of the shape of the output node

    """
    if depth != 0:
        output_image_shape = np.divide(image_shape, np.multiply(pool_size, depth)).tolist()
    else:
        output_image_shape = image_shape
    return tuple([None, filters] + [int(x) for x in output_image_shape])


def get_upconv(depth, nb_filters, pool_size, image_shape, kernel_size=(2, 2, 2), strides=(2, 2, 2),
               deconvolution=False):
    if deconvolution:
        input_shape = compute_level_output_shape(filters=nb_filters, depth=depth + 1, pool_size=pool_size,
                                                 image_shape=image_shape)
        output_shape = compute_level_output_shape(filters=nb_filters, depth=depth, pool_size=pool_size,
                                                  image_shape=image_shape)
        return Deconvolution3D(filters=nb_filters, kernel_size=kernel_size, output_shape=output_shape, strides=strides,
                               input_shape=input_shape)
    else:
        return UpSampling3D(size=pool_size)


#
# def unet_model_2d_noBN(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2), n_labels=0,
#                        loss='mean_squared_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True,
#                        num_gpus=1, num_outputs=1):
#     if n_labels > 0:
#         is_seg_network = True
#     else:
#         is_seg_network = False
#
#     dim = len(input_shape)
#     if dim == 3:
#         ConvL = Conv2D
#         MaxPoolingL = MaxPooling2D
#         pool_size = (2, 2)
#         UpSamplingL = UpSampling2D
#         filter_shape = (5, 5)
#         out_filter_shape = (1, 1)
#
#     elif dim == 4:
#         ConvL = Conv3D
#         MaxPoolingL = MaxPooling3D
#         pool_size = (2, 2, 2)
#         UpSamplingL = UpSampling3D
#         filter_shape = (3, 3, 3)
#         out_filter_shape = (1, 1, 1)
#
#     print('out filter shape is ' + str(out_filter_shape))
#     # input_shape_list = list(input_shape)
#     # input_shape_list.append(1)
#     # input_shape_append = tuple(input_shape_list)
#     # print(input_shape_append)
#     input_img = Input(shape=input_shape, name='input')
#     convs = []
#     pools = []
#     inputs = []
#     centered_inputs = []
#     endpoints = []
#
#     print('unet depth is ')
#     print(unet_depth)
#     for i in range(unet_depth):
#
#         prev = input_img if i == 0 else pools[i - 1]
#         print(int(num_filters * (2 ** i) / downsize_filters_factor))
#         conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
#                      activation='relu', padding='same', kernel_initializer="he_normal",
#                      name=('conv3D_D_1_%d' % (i)))(prev)
#         # conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
#         conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
#                      activation='relu', padding='same', kernel_initializer="he_normal",
#                      name=('conv3D_D_2_%d' % (i)))(conv)
#         # conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
#         if i < (unet_depth - 1):
#             pools.append(MaxPoolingL(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))
#
#         convs.append(conv)
#
#     for i in range(unet_depth - 1):
#         index = i + unet_depth - 1
#         level = unet_depth - (i + 2)
#         up = concatenate([UpSamplingL(size=pool_size, name=('upsampling_U_%d' % (level + 1)))(convs[index]),
#                           convs[level]], axis=-1, name=('concat_%d' % (level)))
#         conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
#                      kernel_initializer="he_normal",
#                      name=('conv3D_U_1_%d' % (level))
#                      )(up)
#         # conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
#         conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
#                      kernel_initializer="he_normal",
#                      name=('conv3D_U_2_%d' % (level)))(conv)
#         # convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))
#         convs.append(conv)
#
#     # conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
#     # conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu", kernel_initializer="he_normal")(conv)
#     # conv = BatchNormalization()(conv)
#     # center_input = Cropping3D(cropping=(0, 0, 0))(input_img)
#
#     inputs.append(input_img)
#     # centered_inputs.append(center_input)
#     print(convs)
#     endpoints.append(convs[-1])
#
#     up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
#     print(loss)
#     print('is_seg_network' + str(is_seg_network))
#     if is_seg_network == False:
#         print(loss)
#         conv = ConvL(num_outputs, out_filter_shape, activation='linear', name='final_conv_3d')(up)
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
#                 return model, parallel_model
#         else:
#             model = Model(inputs=inputs, outputs=conv)
#             if loss == 'grad_loss':
#                 print(loss)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=grad_loss)
#             elif loss == 'warp_image_loss':
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=warp_image_loss)
#             else:
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss)
#             return model, model
#     else:
#         print('segmentation network')
#         if n_labels > 1:
#             conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d')(up)
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                 return model, model
#         else:
#             conv = Conv3D(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(up)
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                 return model, model


# def class_net(feature_shape, dim, unet_num_filters, n_labels, initial_learning_rate, loss='binary_crossentropy'):
#     if dim == 2:
#         ConvL = Conv2D
#         MaxPoolingL = MaxPooling2D
#         filter_shape = (3, 3)
#         pool_shape = (2, 2)
#     elif dim == 3:
#         ConvL = Conv3D
#         MaxPoolingL = MaxPooling3D
#         filter_shape = (3, 3, 3)
#         pool_shape = (2, 2, 2)
#
#     input_shape_list = list(feature_shape)
#     # input_shape_list.append(1)
#     input_shape_append = tuple(input_shape_list)
#     print(input_shape_append)
#
#     model = Sequential()
#     model.add(ConvL(unet_num_filters, filter_shape, padding='same', activation='relu', input_shape=input_shape_append))
#     model.add(ConvL(unet_num_filters, filter_shape, padding='valid', activation='relu'))
#     model.add(MaxPoolingL(pool_size=pool_shape))
#     model.add(Dropout(0.25))
#
#     model.add(ConvL(2 * unet_num_filters, filter_shape, padding='same', activation='relu'))
#     model.add(ConvL(2 * unet_num_filters, filter_shape, activation='relu'))
#     model.add(MaxPoolingL(pool_size=pool_shape))
#     model.add(Dropout(0.25))
#
#     model.add(ConvL((2 ** 2) * unet_num_filters, filter_shape, padding='same', activation='relu'))
#     model.add(ConvL((2 ** 2) * unet_num_filters, filter_shape, activation='relu'))
#     model.add(MaxPoolingL(pool_size=pool_shape))
#     model.add(Dropout(0.25))
#
#     model.add(ConvL((2 ** 3) * unet_num_filters, filter_shape, padding='same', activation='relu'))
#     model.add(ConvL((2 ** 3) * unet_num_filters, filter_shape, activation='relu'))
#     model.add(MaxPoolingL(pool_size=pool_shape))
#     model.add(Dropout(0.25))
#
#     model.add(ConvL((2 ** 4) * unet_num_filters, filter_shape, padding='same', activation='relu'))
#     model.add(ConvL((2 ** 4) * unet_num_filters, filter_shape, activation='relu'))
#     model.add(MaxPoolingL(pool_size=pool_shape))
#     model.add(Dropout(0.25))
#
#     model.add(Flatten())
#     model.add(Dense(512, activation='relu'))
#     model.add(Dropout(0.5))
#
#     if n_labels == 2:
#         model.add(Dense(1, activation='sigmoid'))
#     else:
#         model.add(Dense(n_labels, activation='softmax'))
#
#     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=loss, metrics=['accuracy'])
#     return model

# def unet_model_2d(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2), n_labels=0,
#                   loss='mean_squared_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True,
#                   num_gpus=1, num_outputs=1):
#     if n_labels > 0:
#         is_seg_network = True
#     else:
#         is_seg_network = False
#
#     dim = len(input_shape)
#     if dim == 2:
#         ConvL = Conv2D
#         MaxPoolingL = MaxPooling2D
#         pool_size = (2, 2)
#         UpSamplingL = UpSampling2D
#         filter_shape = (5, 5)
#         out_filter_shape = (1, 1)
#     elif dim == 3:
#         ConvL = Conv3D
#         MaxPoolingL = MaxPooling3D
#         pool_size = (2, 2, 2)
#         UpSamplingL = UpSampling3D
#         filter_shape = (3, 3, 3)
#         out_filter_shape = (1, 1, 1)
#
#     print('out filter shape is ' + str(out_filter_shape))
#     input_shape_list = list(input_shape)
#     input_shape_list.append(1)
#     input_shape_append = tuple(input_shape_list)
#     print(input_shape_append)
#     input_img = Input(shape=input_shape_append, name='input')
#     convs = []
#     pools = []
#     inputs = []
#     centered_inputs = []
#     endpoints = []
#
#     print('unet depth is ' + unet_depth)
#     for i in range(unet_depth):
#         prev = input_img if i == 0 else pools[i - 1]
#         print(int(num_filters * (2 ** i) / downsize_filters_factor))
#         conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
#                      activation='relu', padding='same', kernel_initializer="he_normal",
#                      name=('conv3D_D_1_%d' % (i)))(prev)
#         conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
#         conv = ConvL(int(num_filters * (2 ** i) / downsize_filters_factor), filter_shape,
#                      activation='relu', padding='same', kernel_initializer="he_normal",
#                      name=('conv3D_D_2_%d' % (i)))(conv)
#         conv = BatchNormalization(name=('bnorm_D_2_%d' % (i)))(conv)
#         if i < (unet_depth - 1):
#             pools.append(MaxPoolingL(pool_size, name=('pool_D_%d' % (i)), data_format='channels_last')(conv))
#         convs.append(conv)
#
#     for i in range(unet_depth - 1):
#         index = i + unet_depth - 1
#         level = unet_depth - (i + 2)
#         up = concatenate([UpSamplingL(size=pool_size, name=('upsampling_U_%d' % (level + 1)))(convs[index]),
#                           convs[level]], axis=-1, name=('concat_%d' % (level)))
#         conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
#                      kernel_initializer="he_normal", name=('conv3D_U_1_%d' % (level)))(up)
#         conv = BatchNormalization(name=('bnorm_U_1_%d' % (level)))(conv)
#         conv = ConvL(num_filters * (2 ** level), filter_shape, padding="same", activation="relu",
#                      kernel_initializer="he_normal", name=('conv3D_U_2_%d' % (level)))(conv)
#         convs.append(BatchNormalization(name=('bnorm_U_2_%d' % (level)))(conv))
#
#     inputs.append(input_img)
#     # centered_inputs.append(center_input)
#     endpoints.append(convs[-1])
#     up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
#
#     print(convs)
#     print(loss)
#     print('is_seg_network' + str(is_seg_network))
#
#     if is_seg_network == False:
#         print(loss)
#         conv = ConvL(num_outputs, out_filter_shape, activation='linear', name='final_conv_3d')(up)
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
#     else:
#         print('segmentation network')
#         if n_labels > 1:
#             conv = ConvL(n_labels, out_filter_shape, activation='softmax', name='final_conv_3d')(up)
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                 return model, model
#         else:
#             conv = ConvL(1, (1, 1, 1), activation='sigmoid', name='final_conv_3d')(up)
#             if num_gpus > 1:
#                 with tf.device('/cpu:0'):
#                     model = Model(inputs=inputs, outputs=conv)
#                     parallel_model = multi_gpu_model(model, gpus=num_gpus)
#                     parallel_model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                     return model, parallel_model
#             else:
#                 model = Model(inputs=inputs, outputs=conv)
#                 model.compile(optimizer=Adam(lr=initial_learning_rate), loss=dice_coef_loss2, )
#                 return model, model

# def unet_model_3d(input_shape, num_filters, unet_depth, downsize_filters_factor=1, pool_size=(2, 2, 2), n_labels=0,
#                   loss='mean_absolute_error', initial_learning_rate=0.00001, deconvolution=False, use_patches=True,
#                   num_gpus=1):
#     """
#     Builds the 3D UNet Keras model.
#     :param input_shape: Shape of the input data (x_size, y_size, z_size).
#     :param downsize_filters_factor: Factor to which to reduce the number of filters. Making this value larger will
#     reduce the amount of memory the model will need during training.
#     :param pool_size: Pool size for the max pooling operations.
#     :param n_labels: Number of binary labels that the model is learning.
#     :param initial_learning_rate: Initial learning rate for the model. This will be decayed during training.
#     :param deconvolution: If set to True, will use transpose convolution(deconvolution) instead of upsamping. This
#     increases the amount memory required during training.
#     :return: Untrained 3D UNet Model
#     """
#     # channels last, make feature shape from (32,32,32) - > (32,32,32,1)
#
#     if n_labels > 0:
#         is_seg_network = True
#     else:
#         is_seg_network = False
#
#     # input_shape_list = list(input_shape)
#     # input_shape_list.append(1)
#     # input_shape_append = tuple(input_shape_list)
#     print(input_shape)
#     input_img = Input(shape=input_shape, name='input')
#     convs = []
#     pools = []
#     inputs = []
#     centered_inputs = []
#     endpoints = []
#     print('unet depth is ')
#     print(unet_depth)
#     for i in range(unet_depth):
#         prev = input_img if i == 0 else pools[i - 1]
#         print(int(num_filters * (2 ** i) / downsize_filters_factor))
#         conv = Conv3D(int(num_filters * (2 ** i) / downsize_filters_factor), (3, 3, 3),
#                       activation='relu', padding='same', kernel_initializer="he_normal",
#                       name=('conv3D_D_1_%d' % (i)))(prev)
#         conv = BatchNormalization(name=('bnorm_D_1_%d' % (i)))(conv)
#         conv = Conv3D(int(num_filters * (2 ** i) / downsize_filters_factor), (3, 3, 3),
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
#         up = concatenate([UpSampling3D(size=pool_size, name=('upsampling_U_%d' % (level + 1)))(convs[index]),
#                           convs[level]], axis=-1, name=('concat_%d' % (level)))
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
#     # conv = ZeroPadding3D(padding=(1, 1, 1))(convs[-1])
#     #    conv = Conv3D(num_filters * 2, (3, 3, 3), padding="valid", activation="relu",
#     #                  kernel_initializer="he_normal")(conv)
#     #    conv = BatchNormalization()(conv)
#     #    center_input = Cropping3D(cropping=(0, 0, 0))(input_img)
#
#     inputs.append(input_img)
#     #    centered_inputs.append(center_input)
#     print(convs)
#     endpoints.append(convs[-1])
#
#     up = concatenate(inputs + endpoints, axis=-1, name='final_concat')
#     print(loss)
#     print('is_seg_network' + str(is_seg_network))
#     if is_seg_network == False:
#         print(loss)
#         conv = Conv3D(1, (1, 1, 1), activation='relu', name='final_conv_3d')(up)
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
#
#
