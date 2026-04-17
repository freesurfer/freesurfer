from time import time

import numpy as np
import torch
from batchgenerators.transforms.color_transforms import BrightnessMultiplicativeTransform, \
    ContrastAugmentationTransform, GammaTransform
from batchgenerators.transforms.noise_transforms import GaussianNoiseTransform, GaussianBlurTransform
from batchgenerators.transforms.resample_transforms import SimulateLowResolutionTransform
from batchgenerators.transforms.spatial_transforms import SpatialTransform, MirrorTransform
from batchgenerators.transforms.utility_transforms import RemoveLabelTransform, NumpyToTensor
from nnunetv2.training.data_augmentation.custom_transforms.deep_supervision_donwsampling import \
    DownsampleSegForDSTransform2
from nnunetv2.training.data_augmentation.custom_transforms.masking import MaskTransform
from nnunetv2.training.data_augmentation.custom_transforms.region_based_training import \
    ConvertSegmentationToRegionsTransform
from nnunetv2.training.data_augmentation.custom_transforms.transforms_for_dummy_2d import Convert3DTo2DTransform, \
    Convert2DTo3DTransform

if __name__ == '__main__':
    regions = ((1, 2, 3), (2, 3), (3, ))
    do_dummy_2d_data_aug = False
    patch_size = (128, 128, 128)
    rotation_for_DA = (0, 2*np.pi)
    deep_supervision_scales = ((1, 1, 1), (0.5, 0.5, 0.5), (0.25, 0.25, 0.25))

    tr_transforms = []
    if do_dummy_2d_data_aug:
        ignore_axes = (0,)
        tr_transforms.append(Convert3DTo2DTransform())
        patch_size_spatial = patch_size[1:]
    else:
        patch_size_spatial = patch_size
        ignore_axes = None

    tr_transforms.append(SpatialTransform(
        patch_size_spatial, patch_center_dist_from_border=None,
        do_elastic_deform=False, alpha=(0, 0), sigma=(0, 0),
        do_rotation=True, angle_x=rotation_for_DA, angle_y=rotation_for_DA, angle_z=rotation_for_DA,
        p_rot_per_axis=1,  # todo experiment with this
        do_scale=True, scale=(0.7, 1.4),
        border_mode_data="constant", border_cval_data=0, order_data=3,
        border_mode_seg="constant", border_cval_seg=-1, order_seg=1,
        random_crop=False,  # random cropping is part of our dataloaders
        p_el_per_sample=0, p_scale_per_sample=1, p_rot_per_sample=1,
        independent_scale_for_each_axis=False  # todo experiment with this
    ))

    if do_dummy_2d_data_aug:
        tr_transforms.append(Convert2DTo3DTransform())

    tr_transforms.append(GaussianNoiseTransform(p_per_sample=1, p_per_channel=1))
    tr_transforms.append(GaussianBlurTransform((0.5, 1.), different_sigma_per_channel=True, p_per_sample=1,
                                               p_per_channel=1))
    tr_transforms.append(BrightnessMultiplicativeTransform(multiplier_range=(0.75, 1.25), p_per_sample=1))
    tr_transforms.append(ContrastAugmentationTransform(p_per_sample=1, p_per_channel=1))
    tr_transforms.append(SimulateLowResolutionTransform(zoom_range=(0.5, 1), per_channel=True,
                                                        p_per_channel=1,
                                                        order_downsample=0, order_upsample=3, p_per_sample=1,
                                                        ignore_axes=ignore_axes))
    tr_transforms.append(GammaTransform((0.7, 1.5), True, True, retain_stats=True, p_per_sample=1))
    tr_transforms.append(GammaTransform((0.7, 1.5), False, True, retain_stats=True, p_per_sample=1))

    tr_transforms.append(MirrorTransform((0, 1, 2)))

    tr_transforms.append(MaskTransform([0, 1, 2, 3],
                                       mask_idx_in_seg=0, set_outside_to=0))

    tr_transforms.append(RemoveLabelTransform(-1, 0))

    tr_transforms.append(ConvertSegmentationToRegionsTransform(regions, 'seg', 'seg'))

    tr_transforms.append(DownsampleSegForDSTransform2(deep_supervision_scales, 0, input_key='seg',
                                                      output_key='seg'))

    tr_transforms.append(NumpyToTensor(['data', 'seg'], 'float'))

    compute_times = [[] for i in range(len(tr_transforms))]

    torch.set_num_threads(1)
    for iter in range(50):
        print(iter)
        data_dict = {'data': np.random.uniform(size=(1, 4, 128, 128, 128)),
                     'seg': np.round(4.5 * np.random.uniform(size=(1, 1, 128, 128, 128)) - 1, decimals=0).astype(np.int8)}
        for i, t in enumerate(tr_transforms):
            st = time()
            data_dict = t(**data_dict)
            compute_times[i].append(time() - st)

    for t, ct in zip(tr_transforms, compute_times):
        print(t.__class__.__name__, np.mean(ct))
