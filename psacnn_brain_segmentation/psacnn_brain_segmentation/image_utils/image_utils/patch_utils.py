def show_slices(slices):
    """ Function to display row of image slices """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap="gray", origin="lower")



def extract_patches_dev(in_img_file, patch_size, intensity_threshold=0, step=1, img_fg_extents=None):
    import skimage
    import numpy as np
    import nibabel as nib
    from image_utils import intensity_standardize_utils
    in_img = nib.load(in_img_file)
    in_img_data = in_img.get_data()
    in_img_data_wmp = intensity_standardize_utils.wm_peak_normalize(in_img_data)



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
    print('arg length is ')
    print(len(args))

    import numpy as np
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
    # print newshape
    # print newshape.__class__
    patches = np.reshape(patches, newshape)
    return patches, idx_x, idx_y, idx_z

def build_image_from_patches(in_patches, patch_size, idx_x, idx_y, idx_z, padded_img_size, patch_crop_size):
    ''' patch_crop_size depends on the size of the cnn filter. If [3,3,3] then [1,1,1]'''
    import numpy as np
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
    #
    # subsampled_x = np.arange(0, padded_img_size[0], step_size[0])
    # subsampled_y = np.arange(0, padded_img_size[1], step_size[1])
    # subsampled_z = np.arange(0, padded_img_size[2], step_size[2])
    #
    # sub_idx_x, sub_idx_y, sub_idx_z = np.meshgrid(subsampled_x,
    #                                                  subsampled_y,
    #                                                  subsampled_z,
    #                                                  sparse=False, indexing='ij')
    # sub_hash_idxs = sub_idx_x.flatten()*1000000 + sub_idx_y.flatten()*1000 + sub_idx_z.flatten()
    #
    # # now only consider sub_has_idxs that are present in hash_idxs
    # fg_sub_hash_idxs = np.intersect1d(hash_idxs, sub_hash_idxs, assume_unique=True)
    # bool_idxs = np.in1d(hash_idxs, sub_hash_idxs)
    # np.arange(hash_idxs.shape[0])[np.in1d(hash_idxs, sub_hash_idxs)]
    # sub_patch_indices = np.arange(hash_idxs.shape[0])[np.in1d(hash_idxs, sub_hash_idxs)]
    #
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




#
# def build_image_from_patches(in_patches, img_size, patch_size, step, img_fg_extents):
#     import numpy as np
#     out_img_data = np.zeros(img_size)
#     min_I = img_fg_extents[0, 0]
#     max_I = img_fg_extents[0, 1]
#     min_J = img_fg_extents[1, 0]
#     max_J = img_fg_extents[1, 1]
#     min_K = img_fg_extents[2, 0]
#     max_K = img_fg_extents[2, 1]
#     num_rows = in_patches.shape[0]
#     num_columns = in_patches.shape[1]
#
#     in_patches = np.reshape(in_patches, [num_rows, ])
#     out_img_data_fill = out_img_data[min_I:max_I, min_J:max_J, min_K:max_K]
#     patch_size_x = patch_size[0]
#     patch_size_y = patch_size[1]
#     patch_size_z = patch_size[2]
#
#
#     last_patch_index_x = out_img_data_fill.shape[0] - patch_size_x + 1
#     last_patch_index_y = out_img_data_fill.shape[1] - patch_size_y + 1
#     last_patch_index_z = out_img_data_fill.shape[2] - patch_size_z + 1
#
#
#
#     for iter_z in range(patch_size_z):
#         for iter_y in range(patch_size_y):
#             for iter_x in range(patch_size_x):
#                 column_idx =
#                 out_img_data_fill[iter_x:step:last_patch_index_x,iter_y,iter_z] = np.reshape(in_patches
#
#

