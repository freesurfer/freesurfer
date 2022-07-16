import os
import numpy as np
from scipy.ndimage.interpolation import affine_transform

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()


class VAE:
    def __init__(self, atlasDir, transform, imageSize, seed):

        self.sess = tf.Session()
        self.imageSize = imageSize
        VAEModelPath = os.path.join(atlasDir, 'VAE')
        VAEModelFileName = os.path.join(VAEModelPath, 'model.ckpt.meta')
        saver = tf.train.import_meta_graph(VAEModelFileName)
        saver.restore(self.sess, tf.train.latest_checkpoint(VAEModelPath))
        print('VAE lesion model loaded')

        # Get some info about the VAE
        VAEInfo = np.load(os.path.join(atlasDir, 'VAE/VAEInfo.npz'))
        self.net_shape = VAEInfo['net_shape']
        # transformation matrix mapping train VAE coordinates to template coordinates
        trainToTemplateMat = VAEInfo['trainToTemplateMat']

        # Combination of transformation matrixes in order to obtain a subject to VAE train space transformation
        # First from subject space to template space, then from template space to VAE train space
        # When combining transformations the order of the transformations is from right to left.
        self.trainToSubjectMat = transform.as_numpy_array @ trainToTemplateMat
        self.subjectToTrainMat = np.linalg.inv(self.trainToSubjectMat)

        # Set tf seed 
        self.seed = seed

        # Create tf placeholder
        self.lesionPlaceholder = tf.placeholder(tf.float32, [1, self.net_shape[0], self.net_shape[1], self.net_shape[2], 1])
        self.net = self.run_net(self.lesionPlaceholder, self.net_shape)

    # Encoder network
    def get_encoder(self, data):
        epsilon = 1e-3

        graph = tf.get_default_graph()

        # First convolution
        conv_1_bias = graph.get_tensor_by_name('aconv_1/bias:0')
        conv_1_kernel = graph.get_tensor_by_name('aconv_1/kernel:0')
        conv_1 = tf.nn.conv3d(data, conv_1_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_1_bias
        conv_1 = tf.nn.relu(conv_1)
        mu_1 = graph.get_tensor_by_name('abatch_n_e_1/moving_mean:0')
        var_1 = graph.get_tensor_by_name('abatch_n_e_1/moving_variance:0')
        beta_1 = graph.get_tensor_by_name('abatch_n_e_1/beta:0')
        gamma_1 = graph.get_tensor_by_name('abatch_n_e_1/gamma:0')
        batch_n_e_1 = gamma_1 * ((conv_1 - mu_1) / tf.sqrt(var_1 + epsilon)) + beta_1

        # Second convolution
        conv_2_bias = graph.get_tensor_by_name('aconv_2/bias:0')
        conv_2_kernel = graph.get_tensor_by_name('aconv_2/kernel:0')
        conv_2 = tf.nn.conv3d(batch_n_e_1, conv_2_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_2_bias
        conv_2 = tf.nn.relu(conv_2)
        mu_2 = graph.get_tensor_by_name('abatch_n_e_2/moving_mean:0')
        var_2 = graph.get_tensor_by_name('abatch_n_e_2/moving_variance:0')
        beta_2 = graph.get_tensor_by_name('abatch_n_e_2/beta:0')
        gamma_2 = graph.get_tensor_by_name('abatch_n_e_2/gamma:0')
        batch_n_e_2 = gamma_2 * ((conv_2 - mu_2) / tf.sqrt(var_2 + epsilon)) + beta_2

        # Third convolution
        conv_3_bias = graph.get_tensor_by_name('aconv_3/bias:0')
        conv_3_kernel = graph.get_tensor_by_name('aconv_3/kernel:0')
        conv_3 = tf.nn.conv3d(batch_n_e_2, conv_3_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_3_bias
        conv_3 = tf.nn.relu(conv_3)
        mu_3 = graph.get_tensor_by_name('abatch_n_e_3/moving_mean:0')
        var_3 = graph.get_tensor_by_name('abatch_n_e_3/moving_variance:0')
        beta_3 = graph.get_tensor_by_name('abatch_n_e_3/beta:0')
        gamma_3 = graph.get_tensor_by_name('abatch_n_e_3/gamma:0')
        batch_n_e_3 = gamma_3 * ((conv_3 - mu_3) / tf.sqrt(var_3 + epsilon)) + beta_3

        # Fourth convolution
        conv_4_bias = graph.get_tensor_by_name('aconv_4/bias:0')
        conv_4_kernel = graph.get_tensor_by_name('aconv_4/kernel:0')
        conv_4 = tf.nn.conv3d(batch_n_e_3, conv_4_kernel, strides=[1, 2, 2, 2, 1], padding='VALID') + conv_4_bias
        var_4 = graph.get_tensor_by_name('abatch_n_e_4/moving_variance:0')
        beta_4 = graph.get_tensor_by_name('abatch_n_e_4/beta:0')
        gamma_4 = graph.get_tensor_by_name('abatch_n_e_4/gamma:0')
        conv_4 = tf.nn.relu(conv_4)
        mu_4 = graph.get_tensor_by_name('abatch_n_e_4/moving_mean:0')
        batch_n_e_4 = gamma_4 * ((conv_4 - mu_4) / tf.sqrt(var_4 + epsilon)) + beta_4

        # Fifth convolution
        conv_5_bias = graph.get_tensor_by_name('aconv_5/bias:0')
        conv_5_kernel = graph.get_tensor_by_name('aconv_5/kernel:0')
        conv_5 = tf.nn.conv3d(batch_n_e_4, conv_5_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + conv_5_bias
        conv_5 = tf.nn.relu(conv_5)
        mu_5 = graph.get_tensor_by_name('abatch_n_e_5/moving_mean:0')
        var_5 = graph.get_tensor_by_name('abatch_n_e_5/moving_variance:0')
        beta_5 = graph.get_tensor_by_name('abatch_n_e_5/beta:0')
        gamma_5 = graph.get_tensor_by_name('abatch_n_e_5/gamma:0')
        batch_n_e_5 = gamma_5 * ((conv_5 - mu_5) / tf.sqrt(var_5 + epsilon)) + beta_5

        # Convolve to mu e sigma instead of using fully connected
        mu_bias = graph.get_tensor_by_name('amu/bias:0')
        mu_kernel = graph.get_tensor_by_name('amu/kernel:0')
        mu = tf.nn.conv3d(batch_n_e_5, mu_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + mu_bias

        sigma_bias = graph.get_tensor_by_name('asigma/bias:0')
        sigma_kernel = graph.get_tensor_by_name('asigma/kernel:0')
        sigma = tf.nn.conv3d(batch_n_e_5, sigma_kernel, strides=[1, 1, 1, 1, 1], padding='VALID') + sigma_bias
        sigma = tf.nn.softplus(sigma)

        return mu, sigma

    # decoder network
    def get_decoder(self, code, imageSize):
        epsilon = 1e-3

        graph = tf.get_default_graph()

        code_size = tf.shape(code)

        # First deconv layer
        deconv_0_bias = graph.get_tensor_by_name('adeconv_0/bias:0')
        deconv_0_kernel = graph.get_tensor_by_name('adeconv_0/kernel:0')
        # Deconvolution shape for VALID = stride * (input - 1) + kernel size
        deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3,
                                 code_size[3] - 1 + 3, 16])

        hidden_0_dec = tf.nn.conv3d_transpose(code, deconv_0_kernel, output_shape=deconv_shape,
                                              strides=[1, 1, 1, 1, 1], padding='VALID') + deconv_0_bias

        mu_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_mean:0')
        var_1 = graph.get_tensor_by_name('abatch_n_d_0/moving_variance:0')
        beta_1 = graph.get_tensor_by_name('abatch_n_d_0/beta:0')
        gamma_1 = graph.get_tensor_by_name('abatch_n_d_0/gamma:0')
        batch_n_d_1 = gamma_1 * ((hidden_0_dec - mu_1) / tf.sqrt(var_1 + epsilon)) + beta_1

        # Second deconv layer
        code_size = tf.shape(batch_n_d_1)
        deconv_1_bias = graph.get_tensor_by_name('adeconv_1/bias:0')
        deconv_1_kernel = graph.get_tensor_by_name('adeconv_1/kernel:0')
        deconv_shape = tf.stack([code_size[0], code_size[1] - 1 + 3, code_size[2] - 1 + 3,
                                 code_size[3] - 1 + 3, 16])

        hidden_2_dec = tf.nn.conv3d_transpose(batch_n_d_1, deconv_1_kernel, output_shape=deconv_shape,
                                              strides=[1, 1, 1, 1, 1], padding='VALID') + deconv_1_bias

        hidden_2_dec = tf.nn.relu(hidden_2_dec)
        mu_2 = graph.get_tensor_by_name('abatch_n_d_1/moving_mean:0')
        var_2 = graph.get_tensor_by_name('abatch_n_d_1/moving_variance:0')
        beta_2 = graph.get_tensor_by_name('abatch_n_d_1/beta:0')
        gamma_2 = graph.get_tensor_by_name('abatch_n_d_1/gamma:0')
        batch_n_d_2 = gamma_2 * ((hidden_2_dec - mu_2) / tf.sqrt(var_2 + epsilon)) + beta_2

        # Third deconv layer
        code_size = tf.shape(batch_n_d_2)
        deconv_2_bias = graph.get_tensor_by_name('adeconv_2/bias:0')
        deconv_2_kernel = graph.get_tensor_by_name('adeconv_2/kernel:0')

        deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                                 2 * (code_size[3] - 1) + 5, 16])

        hidden_3_dec = tf.nn.conv3d_transpose(batch_n_d_2, deconv_2_kernel, output_shape=deconv_shape,
                                              strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_2_bias

        hidden_3_dec = tf.nn.relu(hidden_3_dec)
        mu_3 = graph.get_tensor_by_name('abatch_n_d_2/moving_mean:0')
        var_3 = graph.get_tensor_by_name('abatch_n_d_2/moving_variance:0')
        beta_3 = graph.get_tensor_by_name('abatch_n_d_2/beta:0')
        gamma_3 = graph.get_tensor_by_name('abatch_n_d_2/gamma:0')
        batch_n_d_3 = gamma_3 * ((hidden_3_dec - mu_3) / tf.sqrt(var_3 + epsilon)) + beta_3

        # Fourth deconv layer
        code_size = tf.shape(batch_n_d_3)
        deconv_3_bias = graph.get_tensor_by_name('adeconv_3/bias:0')
        deconv_3_kernel = graph.get_tensor_by_name('adeconv_3/kernel:0')

        deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                                 2 * (code_size[3] - 1) + 5, 24])

        hidden_4_dec = tf.nn.conv3d_transpose(batch_n_d_3, deconv_3_kernel, output_shape=deconv_shape,
                                              strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_3_bias

        hidden_4_dec = tf.nn.relu(hidden_4_dec)
        mu_4 = graph.get_tensor_by_name('abatch_n_d_3/moving_mean:0')
        var_4 = graph.get_tensor_by_name('abatch_n_d_3/moving_variance:0')
        beta_4 = graph.get_tensor_by_name('abatch_n_d_3/beta:0')
        gamma_4 = graph.get_tensor_by_name('abatch_n_d_3/gamma:0')
        batch_n_d_4 = gamma_4 * ((hidden_4_dec - mu_4) / tf.sqrt(var_4 + epsilon)) + beta_4

        # Fifth deconv layer
        code_size = tf.shape(batch_n_d_4)
        deconv_4_bias = graph.get_tensor_by_name('adeconv_4/bias:0')
        deconv_4_kernel = graph.get_tensor_by_name('adeconv_4/kernel:0')

        deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                                 2 * (code_size[3] - 1) + 5, 32])

        hidden_5_dec = tf.nn.conv3d_transpose(batch_n_d_4, deconv_4_kernel, output_shape=deconv_shape,
                                              strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_4_bias

        hidden_5_dec = tf.nn.relu(hidden_5_dec)
        mu_5 = graph.get_tensor_by_name('abatch_n_d_4/moving_mean:0')
        var_5 = graph.get_tensor_by_name('abatch_n_d_4/moving_variance:0')
        beta_5 = graph.get_tensor_by_name('abatch_n_d_4/beta:0')
        gamma_5 = graph.get_tensor_by_name('abatch_n_d_4/gamma:0')
        batch_n_d_5 = gamma_5 * ((hidden_5_dec - mu_5) / tf.sqrt(var_5 + epsilon)) + beta_5

        # Sixth deconv layer
        code_size = tf.shape(batch_n_d_5)
        deconv_5_bias = graph.get_tensor_by_name('adeconv_5/bias:0')
        deconv_5_kernel = graph.get_tensor_by_name('adeconv_5/kernel:0')
        deconv_shape = tf.stack([code_size[0], 2 * (code_size[1] - 1) + 5, 2 * (code_size[2] - 1) + 5,
                                 2 * (code_size[3] - 1) + 5, 1])

        hidden_6_dec = tf.nn.conv3d_transpose(batch_n_d_5, deconv_5_kernel, output_shape=deconv_shape,
                                              strides=[1, 2, 2, 2, 1], padding='VALID') + deconv_5_bias

        # We put -100 to the padding to be sure that the final prob is almost if not zero
        hidden_6_dec = self.pad_up_to(hidden_6_dec, [1, imageSize[0], imageSize[1], imageSize[2], 1],
                                 constant_values=-100)

        return tf.nn.sigmoid(hidden_6_dec)

    # add paddings when the size last layer does not match the input size, this happens for deconvolution layers
    def pad_up_to(self, t, max_in_dims, constant_values):
        s = tf.shape(t)
        paddings = [[tf.cast((m - s[i]) / 2, tf.int32), (m - s[i]) - tf.cast((m - s[i]) / 2, tf.int32)]
                    for (i, m) in enumerate(max_in_dims)]
        return tf.pad(t, tf.cast(tf.stack(paddings), tf.int32), 'CONSTANT', constant_values=constant_values)

    # here is the net function. The input goes through the encoder, we sample from it and then it goes through the decoder
    def run_net(self, lesion, imageSize):
        mu, sigma = self.get_encoder(lesion)
        sample_latent = tf.random.normal(mu.shape, 0, 1, seed=self.seed) * sigma + mu
        return self.get_decoder(sample_latent, imageSize)

    def sample(self, lesion):
        # We first go from subject space to train space of the VAE
        # Since we are using scipy affine transform that takes an INVERSE transformation
        # we pass to the function the inverse of subjectToTrainMat, so trainToSubjectMat
        lesionTrainSpace = affine_transform(lesion, self.trainToSubjectMat, output_shape=self.net_shape, order=1)

        # We go through the VAE to get the factorized prior
        # tensorflow wants a 5 dimensional input where the first dimension is batch number
        # and last dimension is number of channels, both set to 1 in our case.
        expandedLesion = np.expand_dims(np.expand_dims(lesionTrainSpace, 0), 4)
        lesionVAETrainSpace = np.squeeze(np.squeeze(self.sess.run(self.net, {self.lesionPlaceholder: expandedLesion}), 4), 0)

        # We then go back to subject space from train space
        # Also here, since we are using scipy affine transform that takes an INVERSE transformation
        # we pass to the function the inverse of trainToSubjectMat, so subjectToTrainMat
        lesionPriorVAE = affine_transform(lesionVAETrainSpace, self.subjectToTrainMat, output_shape=self.imageSize, order=1)

        return lesionPriorVAE
