from typing import Tuple, List, Union
import torch

from batchgeneratorsv2.transforms.base.basic_transform import SegOnlyTransform
from torch.nn.functional import interpolate


class DownsampleSegForDSTransform(SegOnlyTransform):
    def __init__(self, ds_scales: Union[List, Tuple]):
        super().__init__()
        self.ds_scales = ds_scales

    def _apply_to_segmentation(self, segmentation: torch.Tensor, **params) -> List[torch.Tensor]:
        results = []
        for s in self.ds_scales:
            if not isinstance(s, (tuple, list)):
                s = [s] * (segmentation.ndim - 1)
            else:
                assert len(s) == segmentation.ndim - 1

            if all([i == 1 for i in s]):
                results.append(segmentation)
            else:
                new_shape = [round(i * j) for i, j in zip(segmentation.shape[1:], s)]
                dtype = segmentation.dtype
                # interpolate is not defined for short etc
                results.append(interpolate(segmentation[None].float(), new_shape, mode='nearest-exact')[0].to(dtype))
        return results


if __name__ == '__main__':
    from time import time
    import numpy as np
    import os

    os.environ['OMP_NUM_THREADS'] = '1'
    torch.set_num_threads(1)

    mbt = DownsampleSegForDSTransform((1, 0.5, 0.25))

    times_torch = []
    for _ in range(1):
        data_dict = {'segmentation': torch.round(5 * torch.rand((2, 128, 192, 64)), decimals=0).to(torch.uint8)}
        st = time()
        out = mbt(**data_dict)
        times_torch.append(time() - st)
    print('torch', np.mean(times_torch))

    from nnunetv2.training.data_augmentation.custom_transforms.deep_supervision_donwsampling import \
        DownsampleSegForDSTransform2

    gnt_bg = DownsampleSegForDSTransform2((1, 0.5, 0.25), order=0)
    times_bg = []
    for _ in range(1):
        data_dict = {'seg': np.round(5 * np.random.uniform(size=(1, 2, 128, 192, 64)), decimals=0).astype(np.uint8)}
        st = time()
        out = gnt_bg(**data_dict)
        times_bg.append(time() - st)
    print('bg', np.mean(times_bg))
