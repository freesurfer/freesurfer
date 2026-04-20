from typing import List

from batchgeneratorsv2.transforms.base.basic_transform import BasicTransform


class MaskImageTransform(BasicTransform):
    def __init__(self,
                 apply_to_channels: List[int],
                 channel_idx_in_seg: int = 0,
                 set_outside_to: float = 0,
                 ):
        super().__init__()
        self.apply_to_channels = apply_to_channels
        self.channel_idx_in_seg = channel_idx_in_seg
        self.set_outside_to = set_outside_to

    def apply(self, data_dict, **params):
        mask = data_dict['segmentation'][self.channel_idx_in_seg] < 0
        for c in range(data_dict['image'].shape[0]):
            data_dict['image'][c][mask] = self.set_outside_to
        return data_dict

