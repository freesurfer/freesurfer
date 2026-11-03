from copy import deepcopy
from typing import Tuple, List, Union

import numpy as np
import pandas as pd
import torch
from scipy.ndimage import fourier_gaussian
from torch import Tensor
from torch.nn.functional import grid_sample

from batchgeneratorsv2.helpers.scalar_type import RandomScalar, sample_scalar
from batchgeneratorsv2.transforms.base.basic_transform import BasicTransform
from batchgeneratorsv2.transforms.utils.cropping import crop_tensor


class SpatialTransform(BasicTransform):
    def __init__(self, patch_size: Tuple[int, ...],
                 patch_center_dist_from_border: Union[int, List[int], Tuple[int, ...]],
                 random_crop: bool,
                 p_elastic_deform: float = 0, elastic_deform_scale: RandomScalar = (0, 0.2),
                 elastic_deform_magnitude: RandomScalar = (0, 0.2),
                 p_synchronize_def_scale_across_axes: float = False,
                 p_rotation: float = 0, rotation: RandomScalar = (0, 2 * np.pi),
                 p_scaling: float = 0, scaling: RandomScalar = (0.7, 1.3), p_synchronize_scaling_across_axes: float = False,
                 bg_style_seg_sampling: bool = True,
                 mode_seg: str = 'bilinear'
                 ):
        super().__init__()
        self.patch_size = patch_size
        if not isinstance(patch_center_dist_from_border, (tuple, list)):
            patch_center_dist_from_border = [patch_center_dist_from_border] * len(patch_size)
        self.patch_center_dist_from_border = patch_center_dist_from_border
        self.random_crop = random_crop
        self.p_elastic_deform = p_elastic_deform
        self.elastic_deform_scale = elastic_deform_scale  # sigma for blurring offsets, in % of patch size. Larger values mean coarser deformation
        self.elastic_deform_magnitude = elastic_deform_magnitude  # determines the maximum displacement, measured in % of patch size
        self.p_rotation = p_rotation
        self.rotation = rotation
        self.p_scaling = p_scaling
        self.scaling = scaling  # larger numbers = smaller objects!
        self.p_synchronize_scaling_across_axes = p_synchronize_scaling_across_axes
        self.p_synchronize_def_scale_across_axes = p_synchronize_def_scale_across_axes
        self.bg_style_seg_sampling = bg_style_seg_sampling
        self.mode_seg = mode_seg

    def get_parameters(self, **data_dict) -> dict:
        dim = data_dict['image'].ndim - 1

        do_rotation = np.random.uniform() < self.p_rotation
        do_scale = np.random.uniform() < self.p_scaling
        do_deform = np.random.uniform() < self.p_elastic_deform

        if do_rotation:
            angles = [sample_scalar(self.rotation, image=data_dict['image'], dim=i) for i in range(dim)]
        else:
            angles = [0] * dim
        if do_scale:
            scales = [sample_scalar(self.scaling, image=data_dict['image'], dim=i) for i in
                      range(dim)] if np.random.uniform() < self.p_synchronize_scaling_across_axes else [sample_scalar(
                self.scaling, image=data_dict['image'], dim=None)] * dim
        else:
            scales = [1] * dim
        # affine matrix
        if do_scale or do_rotation:
            if dim == 3:
                affine = create_affine_matrix_3d(angles, scales)
            elif dim == 2:
                affine = create_affine_matrix_2d(angles[0], scales)
            else:
                raise RuntimeError(f'Unsupported dimension: {dim}')
        else:
            affine = None  # this will allow us to detect that we can skip computations

        # elastic deformation. We need to create the displacement field here
        # we use the method from augment_spatial_2 in batchgenerators
        if do_deform:
            grid_scale = [i / j for i, j in zip(data_dict['image'].shape[1:], self.patch_size)]
            deformation_scales = [
                sample_scalar(self.elastic_deform_scale, image=data_dict['image'], dim=i, patch_size=self.patch_size)
                for i in
                range(dim)] if np.random.uniform() < self.p_synchronize_scaling_across_axes else [sample_scalar(
                self.elastic_deform_scale, image=data_dict['image'], dim=None, patch_size=self.patch_size)] * dim
            # sigmas must be in pixels, as this will be applied to the deformation field
            sigmas = [i * j for i, j in zip(deformation_scales, self.patch_size)]
            # the magnitude of the deformation field must adhere to the torch's value range for grid_sample, i.e. [-1. 1] and not pixel coordinates. Do not use sigmas here
            # we need to correct magnitude by grid_scale to account for the fact that the grid will be wrt to the image size but the magnitude should be wrt the patch size. oof.
            magnitude = [
                sample_scalar(self.elastic_deform_magnitude, image=data_dict['image'], patch_size=self.patch_size,
                              dim=i, deformation_scale=deformation_scales[i]) / grid_scale[i] for i in range(dim)]
            # doing it like this for better memory layout for blurring
            offsets = torch.normal(mean=0, std=1, size=(dim, *self.patch_size))

            # all the additional time elastic deform takes is spent here
            for d in range(dim):
                # fft torch, slower
                # for i in range(offsets.ndim - 1):
                #     offsets[d] = blur_dimension(offsets[d][None], sigmas[d], i, force_use_fft=True, truncate=6)[0]

                # fft numpy, this is faster o.O
                tmp = np.fft.fftn(offsets[d].numpy())
                tmp = fourier_gaussian(tmp, sigmas)
                offsets[d] = torch.from_numpy(np.fft.ifftn(tmp).real)

                mx = torch.max(torch.abs(offsets[d]))
                offsets[d] /= (mx / np.clip(magnitude[d], a_min=1e-8, a_max=np.inf))
            offsets = torch.permute(offsets, (1, 2, 3, 0))
        else:
            offsets = None
        # grid center must be in [-1, 1] as required by grid_sample
        shape = data_dict['image'].shape[1:]
        if not self.random_crop:
            center_location_in_pixels = [i / 2 for i in shape]
        else:
            center_location_in_pixels = []
            for d in range(dim):
                mn = self.patch_center_dist_from_border[d]
                mx = shape[d] - self.patch_center_dist_from_border[d]
                if mx < mn:
                    center_location_in_pixels.append(shape[d] / 2)
                else:
                    center_location_in_pixels.append(np.random.uniform(mn, mx))
        return {
            'affine': affine,
            'elastic_offsets': offsets,
            'center_location_in_pixels': center_location_in_pixels
        }

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        if params['affine'] is None and params['elastic_offsets'] is None:
            # No spatial transformation is being done. Round grid_center and crop without having to interpolate.
            # This saves compute.
            # cropping requires the center to be given as integer coordinates
            img = crop_tensor(img, [round(i) for i in params['center_location_in_pixels']], self.patch_size, pad_mode='constant',
                              pad_kwargs={'value': 0})
            return img
        else:
            grid = _create_identity_grid(self.patch_size)

            # the grid must be scaled. The grid is [-1, 1] in image coordinates, but we want it to represent the smaller patch
            grid_scale = torch.Tensor([i / j for i, j in zip(img.shape[1:], self.patch_size)])
            grid /= grid_scale

            # we deform first, then rotate
            if params['elastic_offsets'] is not None:
                grid += params['elastic_offsets']
            if params['affine'] is not None:
                grid = torch.matmul(grid, torch.from_numpy(params['affine']).float())

            # we center the grid around the center_location_in_pixels. We should center the mean of the grid, not the center position
            mn = grid.mean(dim=list(range(img.ndim - 1)))
            new_center = torch.Tensor(
                [(j / (i / 2) - 1) for i, j in zip(img.shape[1:], params['center_location_in_pixels'])])
            grid += - mn + new_center
            return grid_sample(img[None], grid[None], mode='bilinear', padding_mode="zeros", align_corners=False)[0]

    def _apply_to_segmentation(self, segmentation: torch.Tensor, **params) -> torch.Tensor:
        segmentation = segmentation.contiguous()
        if params['affine'] is None and params['elastic_offsets'] is None:
            # No spatial transformation is being done. Round grid_center and crop without having to interpolate.
            # This saves compute.
            # cropping requires the center to be given as integer coordinates
            segmentation = crop_tensor(segmentation, [round(i) for i in params['center_location_in_pixels']], self.patch_size,
                                       pad_mode='constant', pad_kwargs={'value': 0})
            return segmentation
        else:
            grid = _create_identity_grid(self.patch_size)

            # the grid must be scaled. The grid is [-1, 1] in image coordinates, but we want it to represent the smaller patch
            grid_scale = torch.Tensor([i / j for i, j in zip(segmentation.shape[1:], self.patch_size)])
            grid /= grid_scale

            # we deform first, then rotate
            if params['elastic_offsets'] is not None:
                grid += params['elastic_offsets']
            if params['affine'] is not None:
                grid = torch.matmul(grid, torch.from_numpy(params['affine']).float())

            # we center the grid around the center_location_in_pixels. We should center the mean of the grid, not the center coordinate
            mn = grid.mean(dim=list(range(segmentation.ndim - 1)))
            new_center = torch.Tensor(
                [(j / (i / 2) - 1) for i, j in zip(segmentation.shape[1:], params['center_location_in_pixels'])])
            grid += - mn + new_center

            if self.mode_seg == 'nearest':
                result_seg = grid_sample(
                                segmentation[None].float(),
                                grid[None],
                                mode=self.mode_seg,
                                padding_mode="zeros",
                                align_corners=False
                            )[0].to(segmentation.dtype)
            else:
                result_seg = torch.zeros((segmentation.shape[0], *self.patch_size), dtype=segmentation.dtype)
                if self.bg_style_seg_sampling:
                    for c in range(segmentation.shape[0]):
                        labels = torch.from_numpy(np.sort(pd.unique(segmentation[c].numpy().ravel())))
                        # if we only have 2 labels then we can save compute time
                        if len(labels) == 2:
                            out = grid_sample(
                                    ((segmentation[c] == labels[1]).float())[None, None],
                                    grid[None],
                                    mode=self.mode_seg,
                                    padding_mode="zeros",
                                    align_corners=False
                                )[0][0] >= 0.5
                            result_seg[c][out] = labels[1]
                            result_seg[c][~out] = labels[0]
                        else:
                            for i, u in enumerate(labels):
                                result_seg[c][
                                    grid_sample(
                                        ((segmentation[c] == u).float())[None, None],
                                        grid[None],
                                        mode=self.mode_seg,
                                        padding_mode="zeros",
                                        align_corners=False
                                    )[0][0] >= 0.5] = u
                else:
                    for c in range(segmentation.shape[0]):
                        labels = torch.from_numpy(np.sort(pd.unique(segmentation[c].numpy().ravel())))
                        #torch.where(torch.bincount(segmentation.ravel()) > 0)[0].to(segmentation.dtype)
                        tmp = torch.zeros((len(labels), *self.patch_size), dtype=torch.float16)
                        scale_factor = 1000
                        done_mask = torch.zeros(*self.patch_size, dtype=torch.bool)
                        for i, u in enumerate(labels):
                            tmp[i] = grid_sample(((segmentation[c] == u).float() * scale_factor)[None, None], grid[None],
                                                 mode=self.mode_seg, padding_mode="zeros", align_corners=False)[0][0]
                            mask = tmp[i] > (0.7 * scale_factor)
                            result_seg[c][mask] = u
                            done_mask = done_mask | mask
                        if not torch.all(done_mask):
                            result_seg[c][~done_mask] = labels[tmp[:, ~done_mask].argmax(0)]
                        del tmp
            del grid
            return result_seg.contiguous()

    def _apply_to_regr_target(self, regression_target, **params) -> torch.Tensor:
        return self._apply_to_image(regression_target, **params)

    def _apply_to_keypoints(self, keypoints, **params):
        raise NotImplementedError

    def _apply_to_bbox(self, bbox, **params):
        raise NotImplementedError


def create_affine_matrix_3d(rotation_angles, scaling_factors):
    # Rotation matrices for each axis
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(rotation_angles[0]), -np.sin(rotation_angles[0])],
                   [0, np.sin(rotation_angles[0]), np.cos(rotation_angles[0])]])

    Ry = np.array([[np.cos(rotation_angles[1]), 0, np.sin(rotation_angles[1])],
                   [0, 1, 0],
                   [-np.sin(rotation_angles[1]), 0, np.cos(rotation_angles[1])]])

    Rz = np.array([[np.cos(rotation_angles[2]), -np.sin(rotation_angles[2]), 0],
                   [np.sin(rotation_angles[2]), np.cos(rotation_angles[2]), 0],
                   [0, 0, 1]])

    # Scaling matrix
    S = np.diag(scaling_factors)

    # Combine rotation and scaling
    RS = Rz @ Ry @ Rx @ S
    return RS


def create_affine_matrix_2d(rotation_angle, scaling_factors):
    # Rotation matrix
    R = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                  [np.sin(rotation_angle), np.cos(rotation_angle)]])

    # Scaling matrix
    S = np.diag(scaling_factors)

    # Combine rotation and scaling
    RS = R @ S
    return RS


def _create_identity_grid(size: List[int]) -> Tensor:
    space = [torch.linspace((-s + 1) / s, (s - 1) / s, s) for s in size[::-1]]
    grid = torch.meshgrid(space, indexing="ij")
    grid = torch.stack(grid, -1)
    spatial_dims = list(range(len(size)))
    grid = grid.permute((*spatial_dims[::-1], len(size)))
    return grid


# An alternative way of generating the displacement fieldQ
# def displacement_field(data: torch.Tensor):
#     downscaling_global = np.random.uniform() ** 2 * 4 + 2
#     # local downscaling can vary a bit relative to global
#     granularity_scale_local = np.random.uniform(round(max(downscaling_global - 1.5, 2)),
#                                                 round(downscaling_global + 1.5), size=3)
#
#     B, _, D, H, W = data.size()
#     random_field_size = [round(j / i) for i, j in zip(granularity_scale_local, data.shape[2:])]
#     pool_kernel_size = [min(i // 4 * 2 + 1, round(7 / 4 * downscaling_global) // 2 * 2 + 1) for i in
#                         random_field_size]  # must be odd
#     pool_padding = [(i - 1) // 2 for i in pool_kernel_size]
#     aug1 = F.avg_pool3d(
#         F.avg_pool3d(
#             torch.randn((B, 2, *random_field_size), device=data.device),
#             pool_kernel_size, stride=1, padding=pool_padding),
#         pool_kernel_size, stride=1, padding=pool_padding)


if __name__ == '__main__':
    torch.set_num_threads(1)

    shape = (128, 128, 128)
    patch_size = (128, 128, 128)
    labels = 2


    # seg = torch.rand([i // 32 for i in shape]) * labels
    # seg_up = torch.round(torch.nn.functional.interpolate(seg[None, None], size=shape, mode='trilinear')[0],
    #                      decimals=0).to(torch.int16)
    # img = torch.ones((1, *shape))
    # img[tuple([slice(img.shape[0])] + [slice(i // 4, i // 4 * 2) for i in shape])] = 200


    import SimpleITK as sitk
    # img = camera()
    # seg = None
    img = sitk.GetArrayFromImage(sitk.ReadImage('/media/isensee/raw_data/nnUNet_raw/Dataset137_BraTS2021/imagesTr/BraTS2021_00000_0000.nii.gz'))
    seg = sitk.GetArrayFromImage(sitk.ReadImage('/media/isensee/raw_data/nnUNet_raw/Dataset137_BraTS2021/labelsTr/BraTS2021_00000.nii.gz'))

    patch_size = (192, 192, 192)
    sp = SpatialTransform(
        patch_size=(192, 192, 192),
        patch_center_dist_from_border=[i / 2 for i in patch_size],
        random_crop=True,
        p_elastic_deform=0,
        elastic_deform_magnitude=(0.1, 0.1),
        elastic_deform_scale=(0.1, 0.1),
        p_synchronize_def_scale_across_axes=0.5,
        p_rotation=1,
        rotation=(-30 / 360 * np.pi, 30 / 360 * np.pi),
        p_scaling=1,
        scaling=(0.75, 1),
        p_synchronize_scaling_across_axes=0.5,
        bg_style_seg_sampling=True,
        mode_seg='bilinear'
    )

    data_dict = {'image': torch.from_numpy(deepcopy(img[None])).float()}
    if seg is not None:
        data_dict['segmentation'] = torch.from_numpy(deepcopy(seg[None]))
    # out = sp(**data_dict)
    #
    # view_batch(out['image'], out['segmentation'])

    from time import time
    times = []
    for _ in range(10):
        data_dict = {'image': torch.from_numpy(deepcopy(img[None])).float()}
        if seg is not None:
            data_dict['segmentation'] = torch.from_numpy(deepcopy(seg[None]))
        st = time()
        out = sp(**data_dict)
        times.append(time() - st)
    print(np.median(times))
