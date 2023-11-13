"""
tensorflow/keras losses for spheremorph
"""

import keras.backend as K
import numpy as np
import tensorflow as tf

from .utils import jacobian_2d, spherical_sin, to_tensor, unpad_2d_image


def get_finite_diff_to_neighbors(img):
    """
           (d)
            |
    (b) -- (a) -- (c)
            |
           (e)
    df_l: a - b
    df_r: a - c
    df_t: a - d
    df_b: a - e
    """
    df_l = img - tf.roll(img, 1, axis=2)  # right shift by 1 to get the left neighbor (a - b)
    df_r = img - tf.roll(img, -1, axis=2)  # left shift by 1 to get the right neighbor (a - c)
    df_t = img - tf.roll(img, 1, axis=1)  # bottom shift by 1 to get the top neighbor (a - d)
    df_b = img - tf.roll(img, -1, axis=1)  # top shift by 1 to get bottom neighbor (a - e)

    return df_l, df_r, df_t, df_b


class SphereLoss:
    def __init__(self,
                 image_shape,
                 pad_size,
                 image_std=None,
                 mask=None,
                 mask_ft_flag=None,
                 eps=1e-3,
                 int_downsize=1):
        """
        class for spherical losses
        :param image_shape: the size of the original 2D image without padding
        :param pad_size: the size of pad
        :param image_std: the precomputed fixed atlas standard divation without padding
        :param mask: mask without padding indicating the where the loss should be evaluted spatially
        :param mask_ft_flag: flag along feature dimension to indicate which feature should be masked
        :param eps: the small positive number to avoid division of zero
        :param int_downsize: optional downsampling size for integral
                             (size of the smooth loss depends on this factor)
        """

        self.int_ds = int_downsize
        self.H, self.W = image_shape
        self.H_ds = self.H // self.int_ds
        self.W_ds = self.W // self.int_ds
        self.eps = eps

        self.pad_size = pad_size
        if self.int_ds > 1:
            self.pad_size_ds = self.pad_size // self.int_ds
        else:
            self.pad_size_ds = self.pad_size

        if image_std is not None:
            std = to_tensor(image_std)
            std = tf.where(std >= eps, std, eps)
            self.image_std = std
            self.image_var_inv = 1 / K.square(std)
            self.image_std_ndim = image_std.ndim
        else:
            self.image_std = None
            self.image_var_inv = None
            self.image_std_ndim = 0

        if mask is not None:
            self.mask = to_tensor(mask)
            self.mask_ds = self.downsample_2d_image(self.mask, 'nearest')
            self.mask_ft_flag = mask_ft_flag
            if mask_ft_flag is not None:
                self.num_ft = len(mask_ft_flag)
                self.mask_ndim = 3
                mask_3d_list = list()
                mask_3d_ds_list = list()
                for m in range(self.num_ft):
                    if mask_ft_flag[m] == 1:
                        mask_3d_list.append(self.mask)
                        mask_3d_ds_list.append(self.mask_ds)
                    else:
                        mask_3d_list.append(tf.ones((self.H, self.W), tf.float32))
                        mask_3d_ds_list.append(tf.ones((self.H_ds, self.W_ds), tf.float32))
                self.mask_3d = tf.stack(mask_3d_list, axis=2)
                self.mask_3d_ds = tf.stack(mask_3d_ds_list, axis=2)
            else:
                self.num_ft = -1  # unknown as mask will be applied after feature channel reduction
                self.mask_ndim = 2
                self.mask_3d = self.mask_3d_ds = None
        else:
            self.mask = self.mask_ds = self.mask_3d = self.mask_3d_ds = self.mask_ft_flag = None

        self.S = spherical_sin(self.H, self.W, eps)  # the sin image
        self.S_2 = K.square(self.S)  # sin^2
        self.S_inv = 1 / self.S  # inverse of sin
        self.S_inv_2 = 1 / self.S_2  # inverse of sin^2

        # cache all downsampled version based on int_ds
        self.S_ds = self.downsample_2d_image(self.S)
        self.S_2_ds = self.downsample_2d_image(self.S_2)
        self.S_inv_ds = self.downsample_2d_image(self.S_inv)
        self.S_inv_2_ds = self.downsample_2d_image(self.S_inv_2)

    def l2_loss(self, weight=1.0, y_true_zero=False):
        """
        l2 loss for data fidelity, weighted by image sigma and corrected for spherical distortion
        :param weight: weight for the loss, can be a vector for different features or a scalar
        :param y_true_zero: binary flag indicating if y_true is zero, when there are two inputs
        :return: the loss function
        """

        def loss(*args):
            # ==================== data in 4D [batch, H, W, feature] ====================
            # parse inputs
            y_true, y_pred, atlas_var = self.parse_data_loss_inputs(*args, y_true_zero=y_true_zero)

            # unpad images to the original size
            y_true, y_pred, atlas_var = self.unpad_loss_images(y_true, y_pred, atlas_var)

            # compute l2 mean square error and average over batch (4D -> 3D)
            mse = K.mean(K.square(y_true - y_pred), axis=[0])

            # ==================== data in 3D [H, W, feature] ====================
            # apply variance term
            if atlas_var is not None:
                # if atlas_var is given by the model as a 4D tensor with a batch size of 1
                mse = self.apply_spatial_weights(mse, atlas_var[0, ...], is_div=True)
            else:
                # if variance is precomputed different for each feature, i.e. image_std is 3D
                if (self.image_var_inv is not None) and (self.image_std_ndim == 3):
                    mse = self.apply_spatial_weights(mse, self.image_var_inv)

            # apply mask according to the mask flag (either apply mask or not for each feature)
            if (self.mask_3d is not None) and (self.mask_ndim == 3):
                mse = self.apply_spatial_weights(mse, self.mask_3d)

            # apply different weights to different feature channels
            mse = self.apply_feature_weights(mse, weight, is_scaled=True)

            # the scale above ensures that K.mean here computes a weighted sum
            # average over feature (3D -> 2D)
            mse = K.mean(mse, axis=[-1])

            # ==================== data in 2D [H, W] ====================
            # only when the variance is precomputed, rather than provided by the model
            if atlas_var is None:
                # optionally weigh the norm square by the inverse of variance
                # same for all features (image_std is 2D)
                if (self.image_var_inv is not None) and (self.image_std_ndim == 2):
                    mse = self.apply_spatial_weights(mse, self.image_var_inv)

            # weigh the norm square by the spherical distortion correction (sin(phi))
            mse = self.apply_spatial_weights(mse, self.S)

            # apply mask if mask is valid and same for all features
            if (self.mask is not None) and (self.mask_ndim == 2):
                mse = self.apply_spatial_weights(mse, self.mask)

            # average over spatial dimensions (2D -> 0D)
            mse = K.mean(mse)

            # ==================== data in 0D [scalar] ====================
            # optionally apply a scalar weight
            mse = self.apply_scalar_weight(mse, weight)
            return mse

        return loss

    def correlation_loss(self, weight=1.0, is_corr=False, spatial_weights=None,
                         ft_idx=None, y_true_zero=False):
        """
        Pearson correlation loss for data fidelity,
        weighted by image sigma and corrected for spherical distortion
        :param weight: weight for the loss, can be a vector for different features or a scalar
        :param is_corr: if True, return the correlation coefficient rather than the loss (1 - r)
        :param spatial_weights: customized spatial weights for the loss
        :param ft_idx: feature index to compute the loss for
        :param y_true_zero: binary flag indicating if y_true is zero, when there are two inputs
        """

        def loss(*args):
            # ==================== data in 4D [batch, H, W, feature] ====================
            y_true, y_pred, atlas_var = self.parse_data_loss_inputs(*args, y_true_zero=y_true_zero)

            if ft_idx is not None:
                # useful when tracking the correlation for a specific feature (used in metric)
                y_true = tf.slice(y_true, [0, 0, 0, ft_idx], [-1, -1, -1, 1])
                y_pred = tf.slice(y_pred, [0, 0, 0, ft_idx], [-1, -1, -1, 1])

            # un-pad images to the original size
            y_true, y_pred, atlas_var, sp_w = self.unpad_loss_images(y_true, y_pred,
                                                                     atlas_var, spatial_weights)

            # spherical distortion correction (sin(phi)) on 4D tensor
            S_4d = self.S[tf.newaxis, ..., tf.newaxis]
            y_true = self.apply_spatial_weights(y_true, S_4d)
            y_pred = self.apply_spatial_weights(y_pred, S_4d)

            # apply customized spatial weights if available
            if spatial_weights is not None:
                sp_w_4d = sp_w[tf.newaxis, ..., tf.newaxis]
                y_true = self.apply_spatial_weights(y_true, sp_w_4d)
                y_pred = self.apply_spatial_weights(y_pred, sp_w_4d)

            # normalize the each 2D image to have zero mean and unit norm, dimensions kept
            means_y_true = K.mean(y_true, axis=[1, 2], keepdims=True)
            y_true = l2_normalize(y_true - means_y_true, axis=[1, 2], keepdims=True)
            means_y_pred = K.mean(y_pred, axis=[1, 2], keepdims=True)
            y_pred = l2_normalize(y_pred - means_y_pred, axis=[1, 2], keepdims=True)

            # compute the correlation using the sum of the dot product (4D -> 2D)
            r = K.sum(y_true * y_pred, axis=[1, 2], keepdims=False)

            # ==================== data in 2D [batch, feature] ====================
            # it is VERY important to compute the loss before applying the feature weights
            if not is_corr:
                r = 1 - r  # r is loss now

            # apply different weights to different feature channels
            r = self.apply_feature_weights(r, weight, is_scaled=True)

            # the scale above ensures that K.mean here computes a weighted sum
            # average over batch and feature (3D -> 2D)
            r = K.mean(r)

            # ==================== data in 0D [scalar] ====================
            # optionally apply a scalar weight
            r = self.apply_scalar_weight(r, weight)
            return r

        def l2_normalize(x, axis, keepdims=True):
            # override tf.math.l2_normalize function to avoid the epsilon clip
            sum_sq = K.sum(K.square(x), axis, keepdims)
            norm = K.sqrt(sum_sq)
            return tf.math.divide_no_nan(x, norm)

        return loss

    def dice_loss(self, weight=1.0, mode='label', num_labels=None, ROI_weight=1.0, is_dice=False):
        """
        dice loss for data fidelity and corrected for spherical distortion
        corrently variance and mask are not supported
        :param weight: weight for the loss, can be a vector for different features or a scalar
        :param mode: the input mode, either 'label' or 'prob' (one-hot)
        :param num_labels: number of labels, required when mode is 'label'
        :param ROI_weight: weight for different ROIs
        :param is_dice: if True, return the dice coefficient rather than the loss (1 - dice)
        :return: the loss function
        """
        assert mode in ['label', 'prob'], 'mode must be either "label" or "prob"'
        if mode == 'label':
            tf.debugging.assert_integer(num_labels,
                                        'num_labels must be an integer')
            assert num_labels > 0, 'num_labels must be a positive integer'

        def loss(*args):
            # ==================== data in 4D [batch, H, W, feature] ====================
            # parse inputs
            y_true, y_pred, atlas_var = self.parse_data_loss_inputs(*args, y_true_zero=False)

            # unpad images to the original size
            y_true, y_pred, atlas_var = self.unpad_loss_images(y_true, y_pred, atlas_var)

            if mode == 'label':
                # round the input to integers, especially needed when learning atlas
                # because atlas is updated by gradient in a continuous manner (float)
                y_true = K.cast(K.round(y_true), 'int32')
                y_pred = K.cast(K.round(y_pred), 'int32')
                # convert to one-hot encoding, they will be automatically expanded at the 5th dim
                y_true = K.one_hot(y_true, num_labels)
                y_pred = K.one_hot(y_pred, num_labels)
            elif mode == 'prob':
                # if the input is already prob/one-hot at the last dim,
                # then we need to expand the dim at the 4th dim to be consistent with the label mode
                y_pred = K.expand_dims(y_pred, axis=3)
                y_true = K.expand_dims(y_true, axis=3)

            # ================ data in 5D [batch, H, W, feature, prob/one-hot] ================
            # correct spatial distortion and compute dice coefficient
            S_5d = self.S[tf.newaxis, ..., tf.newaxis, tf.newaxis]
            top = 2 * K.sum(S_5d * y_true * y_pred, axis=[1, 2])
            bottom = K.sum(S_5d * K.square(y_true), axis=[1, 2]) + \
                     K.sum(S_5d * K.square(y_pred), axis=[1, 2])

            # ==================== data in 3D [batch, feature, prob/one-hot] ====================
            dice = tf.math.divide_no_nan(top, bottom)

            # it is VERY important to compute the loss before applying any feature weights
            if not is_dice:
                dice = 1 - dice  # dice is loss now

            # apply different weights to different ROIs using apply_feature_weights
            dice = self.apply_feature_weights(dice, ROI_weight, is_scaled=True)
            # the scale above ensures that K.mean here computes a weighted sum
            # average over ROI (3D -> 2D)
            dice = K.mean(dice, axis=[-1])

            # ==================== data in 2D [batch, feature] ====================
            # apply different weights to different feature channels
            dice = self.apply_feature_weights(dice, weight, is_scaled=True)
            # the scale above ensures that K.mean here computes a weighted sum
            # average over feature (2D -> 1D)
            dice = K.mean(dice, axis=[-1])

            # ==================== data in 1D [batch] ====================
            # average over batch (1D -> 0D)
            dice = K.mean(dice)

            # ==================== data in 0D [scalar] ====================
            # optionally apply a scalar weight
            dice = self.apply_scalar_weight(dice, weight)
            return dice

        return loss

    def jacobian_loss(self, is_max=True):
        """
        jacobian determinant loss for area distortion
        :param is_max: if True, control the maximum jacobian determinant, otherwise the mean
        """

        def loss(x):
            # the 5th dim is optional used for different flows
            # expand to 5D for 4D input (only one flow)
            if tf.rank(x) == 4:
                x = K.expand_dims(x, axis=-1)

            # ==================== data in 5D [batch, H, W, 2, num_flow] ====================
            batch_size = tf.shape(x)[0]
            num_flow = tf.shape(x)[-1]

            # this nested function computes the jacobian for each batch and flow
            # as tensorflow does not support assigning values to a tensor
            def first_dim(m):
                def last_dim(n):
                    return jacobian_2d(x[m, :, :, :, n], det=True, is_replace_nan=True)

                return tf.map_fn(last_dim, tf.range(num_flow), fn_output_signature=tf.float32)

            J = tf.map_fn(first_dim, tf.range(batch_size), fn_output_signature=tf.float32)

            # ==================== data in 4D [batch, num_flow, H, W] ====================
            # the order of dim is changed due to double tf.map_fn, switch back
            J = tf.transpose(J, perm=[0, 2, 3, 1])

            # ==================== data in 4D [batch, H, W, num_flow] ====================
            # unpad the image with possibly downsampled flow size
            # unpad is done after jacobian computation to avoid boundary effect
            J = self.unpad_loss_images(J, is_ds=True)

            # convert jacobian determinant to distortion measure by subtracting 1
            # and take abs to penalize both expansion and contraction
            J = K.abs(J - 1)

            # spherical distortion correction
            J = tf.multiply(J, self.S_ds[tf.newaxis, ..., tf.newaxis])

            # max/average over flows
            if is_max:
                J = K.max(J, axis=-1)
            else:
                J = K.mean(J, axis=-1)

            # ==================== data in 3D [batch, H, W] ====================
            # max/average over spatial dimensions
            if is_max:
                J = K.max(J, axis=[1, 2])
            else:
                J = K.mean(J, axis=[1, 2])

            # ==================== data in 1D [batch, ] ====================
            # average over batch
            return K.mean(J)

        return loss

    def ncc_loss(self, weight=1.0, win_size=9, is_corr=False, ft_idx=None):
        """
        Local (over window) normalized cross correlation with distortion correction
        """

        def loss(I, J):
            ndims = len(tf.shape(I)) - 2
            assert ndims in [1, 2, 3], "volumes should be 1 to 3 dimensions. found: %d" % ndims

            if ft_idx is not None:
                I = tf.slice(I, [0, 0, 0, ft_idx], [-1, -1, -1, 1])
                J = tf.slice(J, [0, 0, 0, ft_idx], [-1, -1, -1, 1])

            # un-pad image to the original size
            if self.pad_size > 0:
                I = unpad_2d_image(I, self.pad_size)
                J = unpad_2d_image(J, self.pad_size)

            num_chn = tf.cast(tf.shape(I)[-1], dtype=tf.float32)
            if isinstance(win_size, (list, tuple)):
                win = win_size
            else:
                win = [win_size, ] * ndims

            # compute CC squares
            I2 = I * I
            J2 = J * J
            IJ = I * J

            # use depthwise convolution kernel so we can apply weights for different channels
            conv_fn = getattr(tf.nn, 'depthwise_conv%dd' % ndims)
            sum_filt = tf.ones([*win, num_chn, 1])
            strides = [1] * (ndims + 2)
            padding = 'SAME'

            # compute local sums via convolution
            I_sum = conv_fn(I, sum_filt, strides, padding)
            J_sum = conv_fn(J, sum_filt, strides, padding)
            I2_sum = conv_fn(I2, sum_filt, strides, padding)
            J2_sum = conv_fn(J2, sum_filt, strides, padding)
            IJ_sum = conv_fn(IJ, sum_filt, strides, padding)

            # compute means (num of elements in window does not scale with the num of channels)
            win_elem = np.prod(win)
            u_I = I_sum / win_elem
            u_J = J_sum / win_elem

            cross = IJ_sum - u_J * I_sum - u_I * J_sum + u_I * u_J * win_elem
            I_var = I2_sum - 2 * u_I * I_sum + u_I * u_I * win_elem
            J_var = J2_sum - 2 * u_J * J_sum + u_J * u_J * win_elem
            cc = tf.math.divide_no_nan(cross, tf.sqrt(I_var * J_var + self.eps))

            # spherical distortion correction (sin(phi))
            S_4D = tf.expand_dims(self.S, axis=0)
            S_4D = tf.expand_dims(S_4D, axis=-1)
            cc = tf.multiply(cc, S_4D)

            # optionally apply different weights to different feature channels
            if isinstance(weight, (int, float)):
                w = [float(weight)]
            else:
                w = np.array(weight)
            w_ts = to_tensor(w)
            if len(w) > 1:
                # always normalize weights to have a sum of one
                w_ts /= K.sum(w_ts)
                # in this case, need to scale up by the number of features
                # so that we can correctly compute the "weighted sum"
                # in the next step by using K.mean
                w_ts *= num_chn

                # broadcast to the feature dim (cc is still 4D at this point)
                cc = tf.multiply(cc, w_ts[tf.newaxis, tf.newaxis, tf.newaxis, :])

            # average over feature and batch
            cc = K.mean(cc)

            # optionally apply a scalar weight
            if len(w) == 1:
                cc = w_ts * cc

            if is_corr:
                return cc
            else:
                # cc is the correlation in [-1, 1], so 1 - cc is the loss in [0, 2] to minimize
                return 1.0 - cc

        return loss

    def gradient_loss(self):
        """
        gradient loss for smoothness, weighted by edge strength and corrected for spherical distortion
        """

        def loss(y_pred):
            # ==================== data in 4D [batch, H, W, 2] ====================
            # compute finite difference to all 4 neighbors
            df_l, _, df_t, _ = get_finite_diff_to_neighbors(y_pred)

            # un-pad image to the original (possibly downsampled) size
            df_l, df_t = self.unpad_loss_images(df_l, df_t, is_ds=True)

            # compute the norm square of the gradient (4D -> 3D)
            df_l_sq = K.sum(K.square(df_l), axis=[-1])
            df_t_sq = K.sum(K.square(df_t), axis=[-1])

            # ==================== data in 3D [batch, H, W] ====================
            # some math here: phi is the variable along longitude (vertical)
            # theta is the variable along latitude (horizontal)
            # grad f = df/d\phi + 1/sin(\phi) * df/d\theta
            # |grad f|^2 = |df/d\phi|^2 + 1/sin^2(\phi) |df/d\theta|^2
            # \int |grad f|^2 d\Omega =
            # \int (|df/d\phi|^2 + 1/sin^2(\phi) |df/d\theta|^2) sin(\phi) d\phi d\theta

            # apply 1/sin^2 weights on the neighbors with same latitude (horizontally)
            df_l_sq = tf.multiply(df_l_sq, self.S_inv_2_ds[tf.newaxis, :, :])

            # sum of the norm square and average over batch (3D -> 2D)
            df_sq = K.mean(df_l_sq + df_t_sq, axis=[0])

            # ==================== data in 2D [H, W] ====================
            # weigh the norm square by the spherical distortion correction
            # (sin for both vertical and horizontal differences)
            df_sq = tf.multiply(df_sq, self.S_ds)

            # apply mask if mask is same for all features
            if (self.mask_ds is not None) and (self.mask_ndim == 2):
                df_sq = tf.multiply(df_sq, self.mask_ds)

            return K.mean(df_sq)

        return loss

    def laplacian_loss(self, weight=1.0):
        """
        laplacian loss for smoothness, weighted by edge strength and corrected for spherical distortion
        the distortion correction for square of laplacian has not been correctlyh implemented
        to be done in the future
        """

        raise NotImplementedError

        def loss(y_true, y_pred):
            # compute finite difference to all 4 neighbors
            df_l, df_r, df_t, df_b = get_finite_diff_to_neighbors(y_pred)

            # un-pad image to the original (possibly downsampled) size
            df_l = unpad_2d_image(df_l, self.pad_size_ds)
            df_r = unpad_2d_image(df_r, self.pad_size_ds)
            df_t = unpad_2d_image(df_t, self.pad_size_ds)
            df_b = unpad_2d_image(df_b, self.pad_size_ds)

            # compute second order diff (plus is used due to the way of computing first order diff)
            # (x - left) + (x - right) = (x - left) - (right - x)
            df_lr = df_l + df_r
            df_tb = df_t + df_b

            # compute the norm square of the laplacian (4D -> 3D)
            df_lr_sq = K.sum(K.square(df_lr), axis=[-1])
            df_tb_sq = K.sum(K.square(df_tb), axis=[-1])

            # apply 1/sin^2 weights on the neighbors with same latitude (horizontally)
            df_lr_sq = tf.multiply(df_lr_sq, self.S_inv_2_ds[tf.newaxis, :, :])

            # sum of the norm square and average over batch (3D -> 2D)
            df_sq = K.mean(df_lr_sq + df_tb_sq, axis=[0])

            # weigh the norm square by the spherical distortion correction
            df_sq = tf.multiply(df_sq, self.S_ds)

            # apply mask if mask is same for all features
            if (self.mask_ds is not None) and (self.mask_ndim == 2):
                df_sq = tf.multiply(df_sq, self.mask_ds)

            return weight * K.mean(df_sq)

        return loss

    def cc_zero_flow_loss(self):
        """
        enforce flows inside cc to be zero to prevent atlas from dilating
        y_pred has a size of [batch_size, H, W, 2]
        """

        def loss(y_true, y_pred):
            if self.mask_ds is not None:
                # un-pad image to the original size
                if self.pad_size > 0:
                    y_pred = unpad_2d_image(y_pred, self.pad_size)

                # compute the magnitude square of the flow (4D -> 3D)
                flow_mag_sq = K.sum(K.square(y_pred), axis=[-1])

                # average over batch (3D -> 2D)
                flow_mag_sq = K.mean(flow_mag_sq, axis=[0])

                # take only the cc region. mask_ds is regions other than cc.
                flow_mag_sq_cc = tf.multiply(flow_mag_sq, 1 - self.mask_ds)

                return K.sum(flow_mag_sq_cc)
            else:
                return 0.0

        return loss

    def downsample_2d_image(self, img, method='bilinear'):
        if self.int_ds > 1:
            img = img[tf.newaxis, :, :, tf.newaxis]
            new_shape = [self.H_ds, self.W_ds]
            img_ds = tf.squeeze(tf.image.resize(img, new_shape, method))
        else:
            img_ds = img

        return img_ds

    @staticmethod
    def parse_data_loss_inputs(*args, y_true_zero=False):
        """
        :param args: can be in one of the four forms:
                     1) y_pred
                     2) y_pred, atlas_var
                     3) y_true, y_pred
                     4) y_true, y_pred, atlas_var
                     in case 2), set is_y_true_zero=True
        :return: separated inputs
        """
        if len(args) == 1:
            y_pred = args[0]
            y_true = K.zeros_like(y_pred)
            atlas_var = None
        elif len(args) == 2:
            if y_true_zero:
                y_pred, atlas_var = args
                y_true = K.zeros_like(y_pred)
            else:
                y_true, y_pred = args
                atlas_var = None
        elif len(args) == 3:
            y_true, y_pred, atlas_var = args
        else:
            raise ValueError('too many number of inputs')
        return y_true, y_pred, atlas_var

    def unpad_loss_images(self, *args, is_ds=False):
        if self.pad_size > 0:
            if is_ds:
                ps = self.pad_size_ds
            elif self.int_ds == 1:
                ps = self.pad_size
            out = [unpad_2d_image(x, ps) if x is not None else None for x in args]
            if len(out) == 1:
                out = out[0]
            return out
        else:
            return args

    def apply_spatial_weights(self, x, spw, is_div=False):
        # assume image x has the same or compatible (broadcastable) shape with spw
        if spw is not None:
            if is_div:
                # division is only used when spw is variance, which is always positive
                spw = tf.where(spw >= self.eps, spw, self.eps)
                spw = 1.0 / spw
            return tf.multiply(x, spw)
        else:
            return x

    def apply_feature_weights(self, x, weights, is_scaled=True):
        # assume the number of elements in weights is the same as number of features in x
        w_ts, w_len = self.regulate_weights(weights)
        if w_len > 1:  # only apply weights if it's a vector
            if is_scaled:
                # scale up by the number of features to correctly compute the "weighted sum"
                # in a later step by using K.mean
                w_ts *= tf.cast(tf.shape(x)[-1], dtype=tf.float32)

            # add extra dimensions at the beginning to broadcast to the (last) feature dimension
            # can't use the following code because tf complains about the shape change of w_ts
            # for _ in range(tf.rank(x) - 1):
            #     w_ts = w_ts[tf.newaxis, :]
            if tf.rank(x) == 2:
                w_ts = w_ts[tf.newaxis, :]
            elif tf.rank(x) == 3:
                w_ts = w_ts[tf.newaxis, tf.newaxis, :]
            elif tf.rank(x) == 4:
                w_ts = w_ts[tf.newaxis, tf.newaxis, tf.newaxis, :]

            x = tf.multiply(x, w_ts)
        return x

    def apply_scalar_weight(self, x, weight):
        w_ts, w_len = self.regulate_weights(weight)
        if w_len == 1:  # only apply weight if it's a scalar
            x = tf.multiply(x, w_ts)
        return x

    @staticmethod
    def regulate_weights(weights):
        # weights can be a scalar or a vector
        # convert it to np array to find out its length
        # convert it to tensor to use it in tensorflow
        if isinstance(weights, (int, float)):
            w = [float(weights)]
        else:
            w = weights
        w = np.array(w)
        w_len = len(w)
        w_ts = to_tensor(w)
        return w_ts, w_len


# only used with end point loss layers, where y_pred is supposed to be a scalar, just return it
def dummy_loss(y_true, y_pred):
    return y_pred
