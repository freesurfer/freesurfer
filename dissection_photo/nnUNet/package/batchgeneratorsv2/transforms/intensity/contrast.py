import torch

from batchgeneratorsv2.helpers.scalar_type import RandomScalar, sample_scalar
from batchgeneratorsv2.transforms.base.basic_transform import ImageOnlyTransform
import numpy as np


class BGContrast():
    def __init__(self, contrast_range):
        self.contrast_range = contrast_range

    def sample_contrast(self, *args, **kwargs):
        if callable(self.contrast_range):
            factor = self.contrast_range()
        else:
            if np.random.random() < 0.5 and self.contrast_range[0] < 1:
                factor = np.random.uniform(self.contrast_range[0], 1)
            else:
                factor = np.random.uniform(max(self.contrast_range[0], 1), self.contrast_range[1])
        return factor

    def __call__(self, *args, **kwargs):
        return self.sample_contrast(*args, **kwargs)

class ContrastTransform(ImageOnlyTransform):
    def __init__(self, contrast_range: RandomScalar, preserve_range: bool, synchronize_channels: bool, p_per_channel: float = 1):
        super().__init__()
        self.contrast_range = contrast_range
        self.preserve_range = preserve_range
        self.synchronize_channels = synchronize_channels
        self.p_per_channel = p_per_channel

    def get_parameters(self, **data_dict) -> dict:
        shape = data_dict['image'].shape
        apply_to_channel = torch.where(torch.rand(shape[0]) < self.p_per_channel)[0]
        if self.synchronize_channels:
            multipliers = torch.Tensor([sample_scalar(self.contrast_range, image=data_dict['image'], channel=None)] * len(apply_to_channel))
        else:
            multipliers = torch.Tensor([sample_scalar(self.contrast_range, image=data_dict['image'], channel=c) for c in apply_to_channel])
        return {
            'apply_to_channel': apply_to_channel,
            'multipliers': multipliers
        }

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        if len(params['apply_to_channel']) == 0:
            return img
        # array notation is not faster, let's leave it like this
        for i in range(len(params['apply_to_channel'])):
            c = params['apply_to_channel'][i]
            mean = img[c].mean()
            if self.preserve_range:
                minm = img[c].min()
                maxm = img[c].max()

            # this is faster than having it in one line because this circumvents reallocating memory
            img[c] -= mean
            img[c] *= params['multipliers'][i]
            img[c] += mean

            if self.preserve_range:
                img[c].clamp_(minm, maxm)

        return img


if __name__ == '__main__':
    from time import time
    import os

    os.environ['OMP_NUM_THREADS'] = '1'
    torch.set_num_threads(1)

    mbt = ContrastTransform(BGContrast((0.75, 1.25)).sample_contrast, True, False, p_per_channel=1)

    times_torch = []
    for _ in range(100):
        data_dict = {'image': torch.ones((2, 128, 192, 64))}
        st = time()
        out = mbt(**data_dict)
        times_torch.append(time() - st)
    print('torch', np.mean(times_torch))

    from batchgenerators.transforms.color_transforms import ContrastAugmentationTransform

    gnt_bg = ContrastAugmentationTransform((0.75, 1.25), preserve_range=True, per_channel=True, p_per_channel=1)
    times_bg = []
    for _ in range(100):
        data_dict = {'data': np.ones((1, 2, 128, 192, 64))}
        st = time()
        out = gnt_bg(**data_dict)
        times_bg.append(time() - st)
    print('bg', np.mean(times_bg))
