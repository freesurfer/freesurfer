from typing import Tuple

import torch

from batchgeneratorsv2.helpers.scalar_type import RandomScalar, sample_scalar
from batchgeneratorsv2.transforms.base.basic_transform import ImageOnlyTransform
from torch.nn.functional import interpolate


class SimulateLowResolutionTransform(ImageOnlyTransform):
    def __init__(self, scale: RandomScalar, synchronize_channels: bool, synchronize_axes: bool,
                 ignore_axes: Tuple[int, ...],
                 allowed_channels: Tuple[int, ...] = None, p_per_channel: float = 1):
        super().__init__()
        self.scale = scale
        self.synchronize_channels = synchronize_channels
        self.synchronize_axes = synchronize_axes
        self.ignore_axes = ignore_axes
        self.allowed_channels = allowed_channels
        self.p_per_channel = p_per_channel

        self.upmodes = {
            1: 'linear',
            2: 'bilinear',
            3: 'trilinear'
        }

    def get_parameters(self, **data_dict) -> dict:
        shape = data_dict['image'].shape
        if self.allowed_channels is None:
            apply_to_channel = torch.where(torch.rand(shape[0]) < self.p_per_channel)[0]
        else:
            apply_to_channel = [i for i in self.allowed_channels if torch.rand(1) < self.p_per_channel]
        if self.synchronize_channels:
            if self.synchronize_axes:
                scales = torch.Tensor([[sample_scalar(self.scale, image=data_dict['image'], channel=None, dim=None)] * (len(shape) - 1)] * len(apply_to_channel))
            else:
                scales = torch.Tensor([[sample_scalar(self.scale, image=data_dict['image'], channel=None, dim=d) for d in range(len(shape) - 1)]] * len(apply_to_channel))
        else:
            if self.synchronize_axes:
                scales = torch.Tensor([[sample_scalar(self.scale, image=data_dict['image'], channel=c, dim=None)]  * (len(shape) - 1) for c in apply_to_channel])
            else:
                scales = torch.Tensor([[sample_scalar(self.scale, image=data_dict['image'], channel=c, dim=d) for d in range(len(shape) - 1)] for c in apply_to_channel])
        if len(scales) > 0:
            scales[:, self.ignore_axes] = 1
        return {
            'apply_to_channel': apply_to_channel,
            'scales': scales
        }

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        orig_shape = img.shape[1:]
        # we cannot batch this because the downsampled shaps will be different for each channel
        for c, s in zip(params['apply_to_channel'], params['scales']):
            new_shape = [round(i * j.item()) for i, j in zip(orig_shape, s)]
            downsampled = interpolate(img[c][None, None], new_shape, mode='nearest-exact')
            img[c] = interpolate(downsampled, orig_shape, mode=self.upmodes[img.ndim - 1])[0, 0]
        return img


if __name__ == '__main__':
    from time import time
    import numpy as np
    import os

    os.environ['OMP_NUM_THREADS'] = '1'
    torch.set_num_threads(1)

    mbt = SimulateLowResolutionTransform((0.1, 1.), synchronize_channels=False, synchronize_axes=False, ignore_axes=None, allowed_channels=None, p_per_channel=1)

    times_torch = []
    for _ in range(30):
        data_dict = {'image': torch.ones((3, 128, 192, 64))}
        st = time()
        out = mbt(**data_dict)
        times_torch.append(time() - st)
    print('torch', np.mean(times_torch))

    from batchgenerators.transforms.resample_transforms import SimulateLowResolutionTransform as SLRT

    gnt_bg = SLRT((0.1, 1), True, p_per_channel=1, order_downsample=0, order_upsample=1, p_per_sample=1)
    times_bg = []
    for _ in range(30):
        data_dict = {'data': np.ones((1, 3, 128, 192, 64))}
        st = time()
        out = gnt_bg(**data_dict)
        times_bg.append(time() - st)
    print('bg', np.mean(times_bg))
