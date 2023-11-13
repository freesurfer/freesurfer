"""
tensorflow/keras layers for spheremorph
"""

import keras.backend as K
import keras.layers as KL
import neurite as ne
import numpy as np
import tensorflow as tf
from keras.layers import Layer

from .utils import unpad_2d_image, pad_2d_image_spherically


class LossEndPoint(Layer):
    """
    Keras Layer: the loss end point that applies loss function to input tensors
    """

    def __init__(self, loss_fn=None, name=None, metric_fn=None, metric_name=None):
        self.loss_fn = loss_fn

        if name is None:
            name = 'lep'

        if isinstance(metric_fn, (list, tuple)):
            self.metric_fn = metric_fn
            self.metric_name = metric_name
        else:
            if metric_fn is not None:
                self.metric_fn = [metric_fn]
                self.metric_name = [metric_name]
            else:
                self.metric_fn = None
                self.metric_name = None

        super(LossEndPoint, self).__init__(name=name)

    def call(self, inputs, **kwargs):
        if self.loss_fn is not None:
            loss = self.loss_fn(*inputs)
        else:
            loss = 0

        if self.metric_fn is not None:
            for m, metric_f in enumerate(self.metric_fn):
                self.add_metric(metric_f(*inputs), name=self.metric_name[m])

        return K.mean(loss)

    def compute_output_shape(self, input_shape):
        return ()


class SphericalLocalParamWithInput(ne.layers.LocalParamWithInput):
    """
    LocalParamWithInput layer with spherical unpadding and padding
    """

    def __init__(self, shape, initializer='RandomNormal', mult=1.0, pad_size=0, **kwargs):
        self.pad_size = pad_size
        super(SphericalLocalParamWithInput, self).__init__(shape, initializer, mult, **kwargs)

    def call(self, x):
        xslice = K.batch_flatten(x)[:, 0:1]
        b = xslice * tf.zeros((1,)) + tf.ones((1,))
        img = K.flatten(self.kernel * self.biasmult)[tf.newaxis, ...]
        y = K.reshape(K.dot(b, img), [-1, *self.shape])
        if self.pad_size > 0:
            y = unpad_2d_image(y, self.pad_size)
            y = pad_2d_image_spherically(y, self.pad_size)
            y.set_shape(self.compute_output_shape(x.shape))
        return y


class MeanStream(Layer):
    # this layer takes a input batch of subject, compute the mean with a forgetting factor

    def __init__(self, forgetting_factor=0.99, **kwargs):
        self.forgetting_factor = forgetting_factor
        super(MeanStream, self).__init__(**kwargs)

    def build(self, input_shape):
        self.mean = self.add_weight(name='mean',
                                    shape=input_shape[1:],
                                    initializer='zeros',
                                    trainable=False)
        super(MeanStream, self).build(input_shape)

    def call(self, x, training=None):
        # only update mean in training mode
        if training:
            # compute the variance for the current batch, the result is a 3D tensor without batch dim
            mean_curr = K.mean(x, axis=[0])
            var_new = self.forgetting_factor * self.mean + (1 - self.forgetting_factor) * mean_curr
            self.mean.assign(var_new)

        # insert new axis on batch to be compatible with subject and atlas mean
        return self.mean[tf.newaxis, ...]

    def compute_output_shape(self, input_shape):
        return 1, *input_shape[1:]


class VarianceStream(Layer):
    """
    Keras Layer: compute and accumulate the variance with a forgetting factor
    """

    def __init__(self, forgetting_factor=0.99, **kwargs):
        self.forgetting_factor = forgetting_factor
        super(VarianceStream, self).__init__(**kwargs)

    def build(self, input_shape):
        subject_batch_shape, _ = input_shape
        self.variance = self.add_weight(name='variance',
                                        shape=subject_batch_shape[1:],
                                        initializer='zeros',
                                        trainable=False)
        super(VarianceStream, self).build(input_shape)

    def call(self, x, training=None):
        # only update variance in training mode
        if training:
            subject_batch, atlas_mean = x
            # compute the variance for the current batch, the result is a 3D tensor without batch dim
            var_curr = K.mean(K.square(subject_batch - atlas_mean), axis=[0])
            var_new = self.forgetting_factor * self.variance + (
                    1 - self.forgetting_factor) * var_curr
            self.variance.assign(var_new)

        # insert new axis on batch to be compatible with subject and atlas mean
        return self.variance[tf.newaxis, ...]

    def compute_output_shape(self, input_shape):
        subject_batch_shape, _ = input_shape
        return (1, *subject_batch_shape[1:])


class ConcatWithPositionalEncoding(KL.Layer):
    """
    concatenate a set of positional encoding tensors to an existing one
    """

    def __init__(self, npos, pad_size=0, **kwargs):
        self.npos = npos
        self.pad_size = pad_size
        super(ConcatWithPositionalEncoding, self).__init__(**kwargs)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'npos': self.npos,
        })
        return config

    def build(self, inshape):  # create the pe tensor
        input_shape = inshape.as_list()[1:-1]
        if self.pad_size > 0:
            input_shape = [x - 2 * self.pad_size for x in input_shape]
        self.pe = _pos_encoding2D(input_shape, self.npos, self.pad_size)
        super().build(inshape)

    def call(self, x):  # concat of x and self.pe
        conc_fn = lambda x: tf.concat([x, self.pe], axis=-1)
        out = tf.map_fn(conc_fn, x)
        return out

    def compute_output_shape(self, input_shape, **kwargs):
        return input_shape[0:-1] + (input_shape[-1] + 2 * self.npos,)


def Stack(axis=-1, **kwargs):
    def stack_func(x):
        return K.stack(x, axis=axis)

    return KL.Lambda(stack_func, **kwargs)


# functional interface of LossEndPoint layer
def create_loss_end(x, loss_fn=None, **kwargs):
    if not isinstance(x, (list, tuple)):
        x = [x]
    return LossEndPoint(loss_fn, **kwargs)(x)


@tf.function
def _pos_encoding2D(inshape, npos, pad_size=0):
    pe = np.zeros(tuple(inshape) + (2 * npos,))

    for pno in range(npos):
        # x axis / horizontal
        wave_len = inshape[1] / (2 ** pno)
        res = 2 * np.pi / wave_len
        pex_1d = np.cos(np.arange(inshape[1]) * res)
        pex_2d = np.repeat(pex_1d[np.newaxis, ...], inshape[0], axis=0)
        pe[..., pno] = pex_2d

        # y axis / vertical
        wave_len = inshape[0] / (2 ** pno - 0.5)  # -0.5 to avoid same value on two poles
        res = 2 * np.pi / wave_len
        pey_1d = np.cos(np.arange(inshape[0]) * res)
        pey_2d = np.repeat(pey_1d[..., np.newaxis], inshape[1], axis=1)
        pe[..., npos + pno] = pey_2d

    if pad_size > 0:
        pe = pad_2d_image_spherically(pe, pad_size, input_no_batch_dim=True)

    return tf.convert_to_tensor(pe, tf.float32)
