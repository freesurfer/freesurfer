import numpy as np
import nibabel as nib
import skimage
import tensorflow.keras.backend as K
from tensorflow.keras.losses import mean_squared_error as mse
import tensorflow as tf

from .. import wm_peak_normalize


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
    """
    Todo:
        Can be combined with the function below via grad_loss(pure=True)

    """
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

    return lambda1 * mse(y_true_zgrad, y_pred_zgrad) + \
           lambda1 * mse(y_true_ygrad, y_pred_ygrad) + \
           lambda1 * mse(y_true_xgrad, y_pred_xgrad)


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

    return mse(y_true, y_pred) + \
           lambda1 * mse(y_true_zgrad, y_pred_zgrad) + \
           lambda1 * mse(y_true_ygrad, y_pred_ygrad) + \
           lambda1 * mse(y_true_xgrad, y_pred_xgrad)


def extract_patches_dev(in_img_file, patch_size, intensity_threshold=0, step=1, img_fg_extents=None):
    in_img = nib.load(in_img_file)
    in_img_data = in_img.get_data()
    in_img_data_wmp = wm_peak_normalize(in_img_data)

    if img_fg_extents is None:
        (I, J, K) = np.where(in_img_data_wmp > intensity_threshold)
        # get min and max extents of row, column, plane
        min_I = np.min(I)
        min_J = np.min(J)
        min_K = np.min(K)
        max_I = np.max(I)
        max_J = np.max(J)
        max_K = np.max(K)
    else:
        min_I = img_fg_extents[0, 0]
        max_I = img_fg_extents[0, 1]
        min_J = img_fg_extents[1, 0]
        max_J = img_fg_extents[1, 1]
        min_K = img_fg_extents[2, 0]
        max_K = img_fg_extents[2, 1]

    in_img_data_wmp_crop = in_img_data_wmp[min_I:max_I, min_J:max_J, min_K:max_K]
    arr_out = skimage.util.view_as_windows(in_img_data_wmp_crop, patch_size, step)
    num_rows = arr_out.shape[0] * arr_out.shape[1] * arr_out.shape[2]
    num_columns = patch_size[0] * patch_size[1] * patch_size[2]

    # store the extents (will be needed to extract corresponding patches from the target)
    in_extents = np.zeros((3,2))
    in_extents[0, 0] = min_I
    in_extents[0, 1] = max_I

    in_extents[1, 0] = min_J
    in_extents[1, 1] = max_J

    in_extents[2, 0] = min_K;
    in_extents[2, 1] = max_K;

    in_patches = arr_out.reshape((num_rows, num_columns))
    return in_patches, in_extents


def extract_patches(in_img_data, patch_size, intensity_threshold, step_size, *args):
    # pad the images by patch_size
    print('arg length is ' + len(args))
    in_img_data_pad = np.pad(in_img_data, ((patch_size[0], patch_size[0]), (patch_size[1], patch_size[1]),
                                           (patch_size[2], patch_size[2])), 'constant', constant_values=0)
    padded_img_size = in_img_data_pad.shape

    if len(args) == 3:
        idx_x = args[0]
        idx_y = args[1]
        idx_z = args[2]
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
        curr_patch = in_img_data_pad[idx_x[patch_iter]:idx_x[patch_iter] + patch_size[0],
                     idx_y[patch_iter]:idx_y[patch_iter] + patch_size[1],
                     idx_z[patch_iter]:idx_z[patch_iter] + patch_size[2]]
        patches.append(curr_patch)
    patches = np.asarray(patches)

    # add channel as a dimension for keras
    newshape = list(patches.shape)
    newshape.append(1)
    patches = np.reshape(patches, newshape)
    return patches, idx_x, idx_y, idx_z


def build_image_from_patches(in_patches, patch_size, idx_x, idx_y, idx_z, padded_img_size, patch_crop_size):
    '''patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]'''
    out_img_data = np.zeros(padded_img_size)
    count_img_data = np.zeros(padded_img_size)
    patch_mask = np.zeros(patch_size)
    patch_mask[0 + patch_crop_size[0]: patch_size[0] - patch_crop_size[0],
               0 + patch_crop_size[1]: patch_size[1] - patch_crop_size[1],
               0 + patch_crop_size[2]: patch_size[2] - patch_crop_size[2]] = 1

    # this code will build an image from a collection of patches sampled at idx_x, idx_y, idx_z indices
    # these idxs are continuous on the image grid (i.e. step_size = 1
    # however we may want to paste a patch after a step_size=5 after current patch
    # averaging the overlap. the snippet below will consider all the fg idxs
    # and consider only the ones which are subsampled at step_size

    # hash_idxs = idx_x*1000000 + idx_y*1000 + idx_z

    # subsampled_x = np.arange(0, padded_img_size[0], step_size[0])
    # subsampled_y = np.arange(0, padded_img_size[1], step_size[1])
    # subsampled_z = np.arange(0, padded_img_size[2], step_size[2])

    # sub_idx_x, sub_idx_y, sub_idx_z = np.meshgrid(subsampled_x,
    #                                                  subsampled_y,
    #                                                  subsampled_z,
    #                                                  sparse=False, indexing='ij')
    # sub_hash_idxs = sub_idx_x.flatten()*1000000 + sub_idx_y.flatten()*1000 + sub_idx_z.flatten()

    # # now only consider sub_has_idxs that are present in hash_idxs
    # fg_sub_hash_idxs = np.intersect1d(hash_idxs, sub_hash_idxs, assume_unique=True)
    # bool_idxs = np.in1d(hash_idxs, sub_hash_idxs)
    # np.arange(hash_idxs.shape[0])[np.in1d(hash_idxs, sub_hash_idxs)]
    # sub_patch_indices = np.arange(hash_idxs.shape[0])[np.in1d(hash_idxs, sub_hash_idxs)]

    # fg_idx_z = np.mod(fg_sub_hash_idxs, 1000)
    # fg_idx_x_y = np.floor_divide(fg_sub_hash_idxs, 1000)
    # fg_idx_y = np.mod(fg_idx_x_y, 1000)
    # fg_idx_x = np.floor_divide(fg_idx_x_y, 1000)

    for patch_iter in range(len(idx_x)):
        out_img_data[idx_x[patch_iter]:idx_x[patch_iter] + patch_size[0],
        idx_y[patch_iter]:idx_y[patch_iter] + patch_size[1],
        idx_z[patch_iter]:idx_z[patch_iter] + patch_size[2]] += np.multiply(np.reshape(in_patches[patch_iter, :],
                                                                           in_patches.shape[1:4]), patch_mask)
        count_img_data[idx_x[patch_iter]:idx_x[patch_iter] + patch_size[0],
        idx_y[patch_iter]:idx_y[patch_iter] + patch_size[1],
        idx_z[patch_iter]:idx_z[patch_iter] + patch_size[2]] += patch_mask

    out_img_data = np.divide(out_img_data, count_img_data)
    out_img_data[np.isnan(out_img_data)] = 0
    # remove the padding
    unpadded_img_size = padded_img_size - np.multiply(patch_size, 2)

    out_img_data = out_img_data[patch_size[0]:patch_size[0] + unpadded_img_size[0],
                    patch_size[1]:patch_size[1] + unpadded_img_size[1],
                    patch_size[2]:patch_size[2] + unpadded_img_size[2]]
    count_img_data = count_img_data[patch_size[0]:patch_size[0] + unpadded_img_size[0],
                    patch_size[1]:patch_size[1] + unpadded_img_size[1],
                    patch_size[2]:patch_size[2] + unpadded_img_size[2]]

    return out_img_data, count_img_data


def fetch_training_data_files(fs_dir, subjects_dir, img_input_type, training_subject_idxs):
    """
    Note: assumes a freesurfer directory structure

    """
    src_subj_dir_list = sorted(glob.glob(os.path.join(fs_dir, subjects_dir, "[!fs]*", "mri")))
    training_data_files = list()
    input_subj_dir_list = list()
    for i in training_subject_idxs:
        input_subj_dir_list.append(src_subj_dir_list[i])
    for src_dir in input_subj_dir_list:
        training_data_files.append(os.path.join(src_dir, img_input_type + ".mgz"))
    return training_data_files
