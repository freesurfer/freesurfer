from time import time

import numpy as np
import torch

from batchgeneratorsv2.transforms.intensity.brightness import MultiplicativeBrightnessTransform
from batchgeneratorsv2.transforms.intensity.contrast import BGContrast, ContrastTransform
from batchgeneratorsv2.transforms.intensity.gamma import GammaTransform
from batchgeneratorsv2.transforms.intensity.gaussian_noise import GaussianNoiseTransform
from batchgeneratorsv2.transforms.noise.gaussian_blur import GaussianBlurTransform
from batchgeneratorsv2.transforms.spatial.low_resolution import SimulateLowResolutionTransform
from batchgeneratorsv2.transforms.spatial.mirroring import MirrorTransform
from batchgeneratorsv2.transforms.spatial.spatial import SpatialTransform
from batchgeneratorsv2.transforms.utils.compose import ComposeTransforms
from batchgeneratorsv2.transforms.utils.deep_supervision_downsampling import DownsampleSegForDSTransform
from batchgeneratorsv2.transforms.utils.nnunet_masking import MaskImageTransform
from batchgeneratorsv2.transforms.utils.pseudo2d import Convert2DTo3DTransform, Convert3DTo2DTransform
from batchgeneratorsv2.transforms.utils.random import RandomTransform
from batchgeneratorsv2.transforms.utils.remove_label import RemoveLabelTansform
from batchgeneratorsv2.transforms.utils.seg_to_regions import ConvertSegmentationToRegionsTransform

if __name__ == '__main__':
    regions = ((1, 2, 3), (2, 3), (3, ))
    do_dummy_2d_data_aug = False
    patch_size = (128, 128, 128)
    rotation_for_DA = (0, 2*np.pi)
    deep_supervision_scales = ((1, 1, 1), (0.5, 0.5, 0.5), (0.25, 0.25, 0.25))

    transforms = []
    if do_dummy_2d_data_aug:
        ignore_axes = (0,)
        transforms.append(Convert3DTo2DTransform())
        patch_size_spatial = patch_size[1:]
    else:
        patch_size_spatial = patch_size
        ignore_axes = None
    transforms.append(
        SpatialTransform(
            patch_size_spatial, patch_center_dist_from_border=0, random_crop=False, p_elastic_deform=0,
            p_rotation=1,
            rotation=rotation_for_DA, p_scaling=1, scaling=(0.7, 1.4), p_synchronize_scaling_across_axes=1
        )
    )
    if do_dummy_2d_data_aug:
        transforms.append(Convert2DTo3DTransform())

    transforms.append(
        GaussianNoiseTransform(
            noise_variance=(0, 0.1),
            p_per_channel=1,
            synchronize_channels=True
        )
    )
    transforms.append(
        GaussianBlurTransform(
            blur_sigma=(0.5, 1.),
            synchronize_channels=False,
            synchronize_axes=False,
            p_per_channel=1, benchmark=True
        ))
    transforms.append(
        MultiplicativeBrightnessTransform(
            multiplier_range=BGContrast((0.75, 1.25)),
            synchronize_channels=False,
            p_per_channel=1
        ))
    transforms.append(
        ContrastTransform(
            contrast_range=BGContrast((0.75, 1.25)),
            preserve_range=True,
            synchronize_channels=False,
            p_per_channel=1
        ))
    transforms.append(
        SimulateLowResolutionTransform(
            scale=(0.5, 1),
            synchronize_channels=False,
            synchronize_axes=True,
            ignore_axes=ignore_axes,
            allowed_channels=None,
            p_per_channel=1
        ))
    transforms.append(
        GammaTransform(
            gamma=BGContrast((0.7, 1.5)),
            p_invert_image=1,
            synchronize_channels=False,
            p_per_channel=1,
            p_retain_stats=1
        ))
    transforms.append(
        GammaTransform(
            gamma=BGContrast((0.7, 1.5)),
            p_invert_image=0,
            synchronize_channels=False,
            p_per_channel=1,
            p_retain_stats=1
        ))
    transforms.append(
        MirrorTransform(
            allowed_axes=(0, 1, 2)
        )
    )

    transforms.append(MaskImageTransform(
        apply_to_channels=[0, 1, 2, 3],
        channel_idx_in_seg=0,
        set_outside_to=0,
    ))

    transforms.append(
        RemoveLabelTansform(-1, 0)
    )

    transforms.append(
        ConvertSegmentationToRegionsTransform(
            regions=regions,
            channel_in_seg=0
        )
    )

    transforms.append(DownsampleSegForDSTransform(ds_scales=deep_supervision_scales))

    compute_times = [[] for i in range(len(transforms))]

    with torch.no_grad():
        torch.set_num_threads(1)
        for iter in range(50):
            print(iter)
            data_dict = {'image': torch.rand((4, 128, 128, 128)),
                         'segmentation': torch.round(4.5 * torch.rand((1, 128, 128, 128)) - 1, decimals=0).to(torch.int8)}
            for i, t in enumerate(transforms):
                st = time()
                data_dict = t(**data_dict)
                compute_times[i].append(time() - st)

    for t, ct in zip(transforms, compute_times):
        print(t.__class__.__name__ if not isinstance(t, RandomTransform) else t.transform.__class__.__name__, np.mean(ct))
