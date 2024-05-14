from typing import Union, List, Tuple

import torch

from batchgeneratorsv2.transforms.base.basic_transform import BasicTransform


class MoveSegAsOneHotToDataTransform(BasicTransform):
    def __init__(self, source_channel_idx: int, all_labels: Union[Tuple[int, ...], List[int]],
                 remove_channel_from_source: bool = True):
        """
        Used in nnU-Net to append segmentations from the previous stage to the image as additional input
        Args:
            source_channel_idx:
            all_labels:
            remove_channel_from_source:
        """
        super().__init__()
        self.source_channel_idx = source_channel_idx
        self.all_labels = all_labels
        self.remove_channel_from_source = remove_channel_from_source

    def apply(self, data_dict, **params):
        seg = data_dict['segmentation'][self.source_channel_idx]
        seg_onehot = torch.zeros((len(self.all_labels), *seg.shape), dtype=data_dict['image'].dtype)
        for i, l in enumerate(self.all_labels):
            seg_onehot[i][seg == l] = 1
        data_dict['image'] = torch.cat((data_dict['image'], seg_onehot))
        if self.remove_channel_from_source:
            remaining_channels = [i for i in range(data_dict['segmentation'].shape[0]) if i != self.source_channel_idx]
            data_dict['segmentation'] = data_dict['segmentation'][remaining_channels]
        return data_dict
