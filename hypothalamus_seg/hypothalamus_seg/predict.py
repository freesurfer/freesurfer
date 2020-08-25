# python imports
import os
import csv
import numpy as np
from copy import deepcopy
from scipy.ndimage import label

# third-party imports
from ext.lab2im import utils
from ext.lab2im import edit_volumes
from ext.neuron import models as nrn_models

import logging
logging.getLogger('tensorflow').disabled = True


def predict(path_images,
            path_segmentations,
            path_model='../data/model.h5',
            path_posteriors=None,
            path_volumes=None):

    assert path_model, "A model file is necessary"
    assert path_segmentations or path_posteriors, "output segmentation (or posteriors) is required"

    # prepare output filepaths
    images_to_segment, path_segmentations, path_posteriors, path_volumes = prepare_output_files(path_images,
                                                                                                path_segmentations,
                                                                                                path_posteriors,
                                                                                                path_volumes)

    # get label and classes lists
    label_list = np.arange(11)

    # prepare volume file if needed
    if path_volumes is not None:
        csv_header = [['subject'] + [str(lab) for lab in label_list[1:]] + ['whole_left'] + ['whole_right']]
        with open(path_volumes, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(csv_header)
        csvFile.close()

    # perform segmentation
    net = None
    previous_model_input_shape = None
    for idx, (path_image, path_segmentation, path_posterior) in enumerate(zip(images_to_segment,
                                                                              path_segmentations,
                                                                              path_posteriors)):
        utils.print_loop_info(idx, len(images_to_segment), 10)

        # preprocess image and get information
        try:
            image, aff, h, im_res, n_channels, n_dims, shape, cropping, crop_idx = preprocess_image(path_image)
            model_input_shape = list(image.shape[1:])
        except Exception as e:
            print('\nthe following problem occured when preprocessing image %s :' % path_image)
            print(e)
            print('resuming program execution\n')
            continue

        # prepare net for first image or if input's size has changed
        if (idx == 0) | (previous_model_input_shape != model_input_shape):

            # check for image size compatibility
            if (idx != 0) & (previous_model_input_shape != model_input_shape):
                print('image of different shape as previous ones, redefining network')
            previous_model_input_shape = model_input_shape
            net = build_model(path_model, model_input_shape, len(label_list))

        # predict posteriors
        try:
            prediction_patch = net.predict(image)
        except Exception as e:
            print('\nthe following problem occured when predicting segmentation of image %s :' % path_image)
            print(e)
            print('\nresuming program execution')
            continue

        # get posteriors and segmentation
        try:
            seg, posteriors = postprocess(prediction_patch, cropping, shape, crop_idx, n_dims, label_list, aff)
        except Exception as e:
            print('\nthe following problem occured when postprocessing segmentation %s :' % path_segmentation)
            print(e)
            print('\nresuming program execution')
            continue

        # compute volumes
        try:
            if path_volumes is not None:
                volumes = np.sum(posteriors[..., 1:], axis=tuple(range(0, len(posteriors.shape) - 1)))
                volumes = np.around(volumes * np.prod(im_res), 3)
                row = [os.path.basename(path_image).replace('.nii.gz', '')] + [str(vol) for vol in volumes]
                row += [np.sum(volumes[:int(len(volumes) / 2)]), np.sum(volumes[int(len(volumes) / 2):])]
                with open(path_volumes, 'a') as csvFile:
                    writer = csv.writer(csvFile)
                    writer.writerow(row)
                csvFile.close()
        except Exception as e:
            print('\nthe following problem occured when computing the volumes of segmentation %s :' % path_segmentation)
            print(e)
            print('\nresuming program execution')
            continue

        # write results to disk
        try:
            if path_segmentation is not None:
                utils.save_volume(seg.astype('int'), aff, h, path_segmentation)
            if path_posterior is not None:
                if n_channels > 1:
                    new_shape = list(posteriors.shape)
                    new_shape.insert(-1, 1)
                    new_shape = tuple(new_shape)
                    posteriors = np.reshape(posteriors, new_shape)
                utils.save_volume(posteriors.astype('float'), aff, h, path_posterior)
        except Exception as e:
            print('\nthe following problem occured when saving the results for image %s :' % path_image)
            print(e)
            print('\nresuming program execution')
            continue

    # print output info
    print('\n')
    if path_segmentations[0] is not None:
        print('segmentations saved in: ' + os.path.dirname(path_segmentations[0]))
    if path_posteriors[0] is not None:
        print('posteriors saved in:    ' + os.path.dirname(path_posteriors[0]))
    if path_volumes is not None:
        print('volumes saved in:       ' + path_volumes)


def prepare_output_files(path_images, out_seg, out_posteriors, out_volumes):

    # convert path to absolute paths
    path_images = os.path.abspath(path_images)
    if out_seg is not None:
        out_seg = os.path.abspath(out_seg)
    if out_posteriors is not None:
        out_posteriors = os.path.abspath(out_posteriors)
    if out_volumes is not None:
        out_volumes = os.path.abspath(out_volumes)

    # prepare input/output volumes
    if ('.nii.gz' not in path_images) & ('.nii' not in path_images) & ('.mgz' not in path_images) & \
            ('.npz' not in path_images):
        images_to_segment = utils.list_images_in_folder(path_images)
        assert len(images_to_segment) > 0, "Could not find any training data"
        if out_seg:
            if not os.path.exists(out_seg):
                os.mkdir(out_seg)
            out_seg = [os.path.join(out_seg, os.path.basename(image)).replace('.nii', '_seg.nii') for image in
                       images_to_segment]
            out_seg = [seg_path.replace('.mgz', '_seg.mgz') for seg_path in out_seg]
            out_seg = [seg_path.replace('.npz', '_seg.npz') for seg_path in out_seg]
        else:
            out_seg = [out_seg] * len(images_to_segment)
        if out_posteriors:
            if not os.path.exists(out_posteriors):
                os.mkdir(out_posteriors)
            out_posteriors = [os.path.join(out_posteriors, os.path.basename(image)).replace('.nii',
                              '_posteriors.nii') for image in images_to_segment]
            out_posteriors = [posteriors_path.replace('.mgz', '_posteriors.mgz')
                              for posteriors_path in out_posteriors]
            out_posteriors = [posteriors_path.replace('.npz', '_posteriors.npz')
                              for posteriors_path in out_posteriors]
        else:
            out_posteriors = [out_posteriors] * len(images_to_segment)

    else:
        assert os.path.exists(path_images), "Could not find image to segment"
        images_to_segment = [path_images]
        if out_seg is not None:
            if ('.nii.gz' not in out_seg) & ('.nii' not in out_seg) & ('.mgz' not in out_seg) & ('.npz' not in out_seg):
                if not os.path.exists(out_seg):
                    os.mkdir(out_seg)
                filename = os.path.basename(path_images).replace('.nii', '_seg.nii')
                filename = filename.replace('mgz', '_seg.mgz')
                filename = filename.replace('.npz', '_seg.npz')
                out_seg = os.path.join(out_seg, filename)
            else:
                if not os.path.exists(os.path.dirname(out_seg)):
                    os.mkdir(os.path.dirname(out_seg))
        out_seg = [out_seg]
        if out_posteriors is not None:
            if ('.nii.gz' not in out_posteriors) & ('.nii' not in out_posteriors) & ('.mgz' not in out_posteriors) & \
                    ('.npz' not in out_posteriors):
                if not os.path.exists(out_posteriors):
                    os.mkdir(out_posteriors)
                filename = os.path.basename(path_images).replace('.nii', '_posteriors.nii')
                filename = filename.replace('mgz', '_posteriors.mgz')
                filename = filename.replace('.npz', '_posteriors.npz')
                out_posteriors = os.path.join(out_posteriors, filename)
            else:
                if not os.path.exists(os.path.dirname(out_posteriors)):
                    os.mkdir(os.path.dirname(out_posteriors))
        out_posteriors = [out_posteriors]

    if out_volumes:
        if out_volumes[-4:] != '.csv':
            print('out_volumes provided without csv extension. Adding csv extension to output_volumes.')
            out_volumes += '.csv'
        if not os.path.exists(os.path.dirname(out_volumes)):
            os.mkdir(os.path.dirname(out_volumes))

    return images_to_segment, out_seg, out_posteriors, out_volumes


def preprocess_image(im_path):

    # read image and corresponding info
    n_levels = 3
    im, shape, aff, n_dims, n_channels, header, im_res = utils.get_volume_info(im_path,
                                                                               aff_ref=np.eye(4),
                                                                               return_volume=True)

    # check that patch_shape or im_shape are divisible by 2**n_levels
    if not all([size % (2**n_levels) == 0 for size in shape]):
        crop_shape = [min(utils.find_closest_number_divisible_by_m(size, 2 ** n_levels), 232) for size in shape]
    else:
        if not all([size <= 232 for size in shape]):
            crop_shape = [min(size, 232) for size in shape]
        else:
            crop_shape = None

    # crop image if necessary
    if crop_shape is not None:
        crop_idx = np.round((shape - np.array(crop_shape)) / 2).astype('int')
        crop_idx = np.concatenate((crop_idx, crop_idx + crop_shape), axis=0)
        im = edit_volumes.crop_volume_with_idx(im, crop_idx=crop_idx)
    else:
        crop_idx = None

    # normalise image
    m = np.min(im)
    M = np.max(im)
    if M == m:
        im = np.zeros(im.shape)
    else:
        im = (im - m) / (M - m)

    # add batch and channel axes
    if n_channels > 1:
        im = utils.add_axis(im)
    else:
        im = utils.add_axis(im, -2)

    return im, aff, header, im_res, n_channels, n_dims, shape, crop_shape, crop_idx


def build_model(model_file, input_shape, n_lab):

    # build UNet
    batch_norm_dim = -1
    net = nrn_models.unet(nb_features=24,
                          input_shape=input_shape,
                          nb_levels=3,
                          conv_size=3,
                          nb_labels=n_lab,
                          name='unet',
                          prefix=None,
                          feat_mult=2,
                          pool_size=2,
                          padding='same',
                          dilation_rate_mult=1,
                          activation='elu',
                          use_residuals=False,
                          final_pred_activation='softmax',
                          nb_conv_per_level=2,
                          layer_nb_feats=None,
                          conv_dropout=0,
                          batch_norm=batch_norm_dim,
                          input_model=None)
    net.load_weights(model_file, by_name=True)

    return net


def postprocess(prediction, crop_shape, im_shape, crop, n_dims, labels, aff):

    # get posteriors and segmentation
    post_patch = np.squeeze(prediction)
    seg_patch = post_patch.argmax(-1)

    # keep biggest connected component (use it with smoothing!)
    seg_left = deepcopy(seg_patch)
    seg_left[seg_left > 5] = 0
    components, n_components = label(seg_left, np.ones([n_dims]*n_dims))
    if n_components > 1:
        unique_components = np.unique(components)
        size = 0
        mask = None
        for comp in unique_components[1:]:
            tmp_mask = components == comp
            tmp_size = np.sum(tmp_mask)
            if tmp_size > size:
                size = tmp_size
                mask = tmp_mask
        seg_left[np.logical_not(mask)] = 0

    seg_right = deepcopy(seg_patch)
    seg_right[seg_right < 6] = 0
    components, n_components = label(seg_right, np.ones([n_dims]*n_dims))
    if n_components > 1:
        unique_components = np.unique(components)
        size = 0
        mask = None
        for comp in unique_components[1:]:
            tmp_mask = components == comp
            tmp_size = np.sum(tmp_mask)
            if tmp_size > size:
                size = tmp_size
                mask = tmp_mask
        seg_right[np.logical_not(mask)] = 0

    seg_patch = seg_left | seg_right

    # align prediction back to first orientation
    seg_patch = edit_volumes.align_volume_to_ref(seg_patch, np.eye(4), aff_ref=aff)
    post_patch = edit_volumes.align_volume_to_ref(post_patch, np.eye(4), aff_ref=aff, n_dims=n_dims)

    # paste patches back to matrix of original image size
    if crop_shape is not None:
        seg = np.zeros(shape=im_shape, dtype='int32')
        posteriors = np.zeros(shape=[*im_shape, labels.shape[0]])
        posteriors[..., 0] = np.ones(im_shape)  # place background around patch
        if n_dims == 2:
            seg[crop[0]:crop[2], crop[1]:crop[3]] = seg_patch
            posteriors[crop[0]:crop[2], crop[1]:crop[3], :] = post_patch
        elif n_dims == 3:
            seg[crop[0]:crop[3], crop[1]:crop[4], crop[2]:crop[5]] = seg_patch
            posteriors[crop[0]:crop[3], crop[1]:crop[4], crop[2]:crop[5], :] = post_patch
    else:
        seg = seg_patch
        posteriors = post_patch
    seg = labels[seg.astype('int')].astype('int')

    return seg, posteriors
