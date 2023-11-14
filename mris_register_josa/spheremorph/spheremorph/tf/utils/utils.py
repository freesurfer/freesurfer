"""
tensorflow/keras utilities for spheremorph
"""

import numpy as np
import tensorflow as tf


def pad_2d_image_spherically(img, pad_size=8, input_no_batch_dim=False):
    """
    pad parameterized 2d image based on the spherical positions of its vertices
    img: image to pad, whose shape is [batch_size, H, W, ...] or [H, W] for a single image
    """
    is_2d = is_nd(img, 2)
    img = expand_batch_dims_with_cond(img, is_2d, input_no_batch_dim)

    if pad_size > 0:
        # pad the north pole on top
        top = img[:, 1:pad_size + 1, ...]  # get top pad without the first row (reflect)
        top = flip(top, axis=1)  # flip upside down
        top = roll(top, get_shape(top, 2) // 2, axis=2)  # circularly shift by pi

        # similarly for the south pole on bottom
        bot = img[:, -pad_size - 1:-1, ...]
        bot = flip(bot, axis=1)
        bot = roll(bot, get_shape(bot, 2) // 2, axis=2)

        # concatenate top and bottom before padding left and right
        img2 = concat((top, img, bot), axis=1)

        # pad left to right and right to left (wrap)
        left = img2[:, :, 0:pad_size, ...]
        right = img2[:, :, -pad_size:, ...]
        img3 = concat((right, img2, left), axis=2)
    else:
        img3 = img

    img3 = squeeze_with_cond(img3, is_2d, input_no_batch_dim)

    return img3


def unpad_2d_image(img, pad_size=0, input_no_batch_dim=False):
    """
    extract the original image from the padded image
    img: image to unpad, whose shape is [batch_size, H, W, ...]
    """
    is_2d = is_nd(img, 2)
    img = expand_batch_dims_with_cond(img, is_2d, input_no_batch_dim)

    if pad_size > 0:
        img = img[:, pad_size:-pad_size, pad_size:-pad_size, ...]

    img = squeeze_with_cond(img, is_2d, input_no_batch_dim)
    return img


# wrap conditional functions compatible with both numpy and tf only for padding and unpadding functions
def expand_batch_dims_with_cond(data, is_2d, input_no_batch_dim):
    if tf.is_tensor(data):
        is_2d = tf.cast(is_2d, tf.bool)
        input_no_batch_dim = tf.cast(input_no_batch_dim, tf.bool)
        cond = tf.logical_or(is_2d, input_no_batch_dim)
        data = tf.cond(cond, lambda: expand_batch_dims(data), lambda: data)
    else:
        if is_2d or input_no_batch_dim:
            data = expand_batch_dims(data)
    return data


def squeeze_with_cond(data, is_2d, input_no_batch_dim):
    if tf.is_tensor(data):
        is_2d = tf.cast(is_2d, tf.bool)
        input_no_batch_dim = tf.cast(input_no_batch_dim, tf.bool)
        cond = tf.logical_or(is_2d, input_no_batch_dim)
        data = tf.cond(cond, lambda: squeeze(data), lambda: data)
    else:
        if is_2d or input_no_batch_dim:
            data = squeeze(data)
    return data


# wrap some basic functions so they are compatible with both numpy and tf
def get_shape(data, axis=None):
    if tf.is_tensor(data):
        if axis is None:
            return tf.shape(data)
        else:
            return tf.shape(data)[axis]
    else:
        if axis is None:
            return data.shape
        else:
            return data.shape[axis]


def concat(data, axis):
    if tf.is_tensor(data[0]):
        return tf.concat(data, axis)
    else:
        return np.concatenate(data, axis)


def flip(data, axis):
    if tf.is_tensor(data):
        if not (isinstance(axis, list) or isinstance(axis, tuple)):
            axis = [axis]
        return tf.reverse(data, axis)
    else:
        return np.flip(data, axis)


def roll(data, shift, axis):
    if tf.is_tensor(data):
        return tf.roll(data, shift, axis)
    else:
        return np.roll(data, shift, axis)


def squeeze(data):
    if tf.is_tensor(data):
        return tf.squeeze(data)
    else:
        return np.squeeze(data)


def expand_batch_dims(data):
    if tf.is_tensor(data):
        return tf.expand_dims(data, 0)
    else:
        return data[np.newaxis, ...]


def is_nd(data, n):
    if tf.is_tensor(data):
        return tf.cond(tf.cast(tf.rank(data) == n, tf.bool),
                       lambda: tf.constant(True, dtype=tf.bool),
                       lambda: tf.constant(False, dtype=tf.bool))
    else:
        if data.ndim == n:
            return True
        else:
            return False


def to_tensor(data, d_type=tf.float32):
    if tf.is_tensor(data):
        if data.dtype is not d_type:
            data = tf.cast(data, d_type)
        return data
    else:
        return tf.convert_to_tensor(data, d_type)


def spherical_sin(H, W, eps=1e-3):
    # the sin of latitude assumes the first and last element hits the poles
    # i.e. sin(0) and sin(pi), so that it can be padded at the north
    # and south poles properly by "reflection"
    rg = tf.transpose(tf.range(0, H, dtype=tf.float32))
    rg = tf.math.divide_no_nan(rg, to_tensor(H - 1))
    rg = tf.math.multiply(rg, to_tensor(np.pi))
    sin_lat = tf.math.sin(rg)
    # sin_lat = tf.math.sin(to_tensor(np.arange(0, H).T / (H - 1) * np.pi))
    # remove singularity near the two poles by setting the minimum to (positive) eps
    sin_lat = tf.where(sin_lat >= eps, sin_lat, eps)
    # repeat on longitude, forming the sin matrix
    S = tf.expand_dims(sin_lat, 1) * tf.ones((1, W))
    return S


def jacobian_2d(x, det=False, is_use_tf_det=False, is_replace_nan=False):
    """
    this function compute jacobian or its determinant of a 2D vector field
    all lines are written in an explicit and tf compatible way so it can be used for loss evaluation
    i.e. it can be excuted in graph mode
    when used for loss evalution, do NOT use tf.linalg.det to compute the determinant as it will
    cause matrix "not invertible" error during backpropagation
    this function was adapted from nes.utils.gradient/jacobian and vxm.py.utils.jacobian_determinant
    :param x: input tensor of shape [B, H, W, 2]
    :param det: if True, return the determinant of the Jacobian
    :param is_use_tf_det: if True, use tf.linalg.det to compute the determinant
    :param is_replace_nan: if True, replace nan values with 1
    :return: Jacobian tensor of shape [B, H, W, 2, 2] or [B, H, W] if det=True
    """

    def gradient_2d(x):
        beg = (x[1, ...] - x[0, ...])[tf.newaxis, ...]
        mid = 0.5 * (x[2:, ...] - x[:-2, ...])
        end = (x[-1, ...] - x[-2, ...])[tf.newaxis, ...]
        dx = tf.concat((beg, mid, end), axis=0)

        beg = (x[:, 1, ...] - x[:, 0, ...])[:, tf.newaxis, ...]
        mid = 0.5 * (x[:, 2:, ...] - x[:, :-2, ...])
        end = (x[:, -1, ...] - x[:, -2, ...])[:, tf.newaxis, ...]
        dy = tf.concat((beg, mid, end), axis=1)

        return dx, dy

    shape = tf.shape(x)
    gridx = tf.range(shape[0], dtype=tf.float32)
    gridy = tf.range(shape[1], dtype=tf.float32)
    grid = tf.meshgrid(gridx, gridy, indexing='ij')
    x += tf.stack(grid, axis=-1)

    grad = gradient_2d(x)

    if det:
        if is_use_tf_det:
            # using tf.linalg.det will replace nan/inf with 0
            J = tf.stack(grad, axis=-1)
            J = tf.linalg.det(J)
        else:
            dfdx = grad[0]
            dfdy = grad[1]
            J = dfdx[..., 0] * dfdy[..., 1] - dfdy[..., 0] * dfdx[..., 1]
            if is_replace_nan:
                # we replace nan/inf with 1 if we want to use J for loss
                J = tf.where(tf.math.is_finite(J), J, tf.ones_like(J))
    else:
        J = tf.stack(grad, axis=-1)

    return J
