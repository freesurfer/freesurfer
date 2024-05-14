from time import time
from typing import Union, List, Tuple
import numpy as np
import torch
from acvl_utils.morphology.morphology_helper import label_with_component_sizes

from batchgeneratorsv2.transforms.base.basic_transform import ImageOnlyTransform


class RemoveRandomConnectedComponentFromOneHotEncodingTransform(ImageOnlyTransform):
    def __init__(self,
                 channel_idx: Union[int, List[int], Tuple[int, ...]],
                 fill_with_other_class_p: float = 0.25,
                 dont_do_if_covers_more_than_x_percent: float = 0.25,
                 p_per_label: float = 1
                 ):
        super().__init__()
        if not isinstance(channel_idx, (list, tuple)):
            channel_idx = [channel_idx]
        if isinstance(channel_idx, tuple):
            channel_idx = list(channel_idx)
        self.channel_idx = channel_idx
        self.fill_with_other_class_p = fill_with_other_class_p
        self.dont_do_if_covers_more_than_x_percent = dont_do_if_covers_more_than_x_percent
        self.p_per_label = p_per_label

    def get_parameters(self, **data_dict) -> dict:
        # this needs to be applied in random order to the channels
        np.random.shuffle(self.channel_idx)
        apply_to_channels = [self.channel_idx[i] for i, j in enumerate(torch.rand(len(self.channel_idx)) < self.p_per_label) if j]

        # self.fill_with_other_class_p cannot be resolved here because we don't know how many components there are
        return {
            'apply_to_channels': apply_to_channels,
        }

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        for a in params['apply_to_channels']:
            workon = img[a].to(bool).numpy()
            if not np.any(workon):
                continue
            num_voxels = np.prod(workon.shape, dtype=np.uint64)
            lab, component_sizes = label_with_component_sizes(workon.astype(bool))
            if len(component_sizes) > 0:
                valid_component_ids = [i for i, j in component_sizes.items() if j <
                                       num_voxels * self.dont_do_if_covers_more_than_x_percent]
                # print('RemoveRandomConnectedComponentFromOneHotEncodingTransform', c,
                # np.unique(data[b, c]), len(component_sizes), valid_component_ids,
                # len(valid_component_ids))
                if len(valid_component_ids) > 0:
                    random_component = np.random.choice(valid_component_ids)
                    img[a][lab == random_component] = 0
                    if np.random.uniform() < self.fill_with_other_class_p:
                        other_ch = [i for i in self.channel_idx if i != a]
                        if len(other_ch) > 0:
                            other_class = np.random.choice(other_ch)
                            img[other_class][lab == random_component] = 1
        return img


if __name__ == '__main__':
    torch.set_num_threads(1)

    tr = RemoveRandomConnectedComponentFromOneHotEncodingTransform(
        channel_idx=(0, 1),
        fill_with_other_class_p=1,
        dont_do_if_covers_more_than_x_percent=0.25,
        p_per_label=1
    )

    times_torch = []
    for _ in range(10):
        # img = (torch.rand((1, 128, 128, 128)) < 0.5).to(torch.int16)
        # img = torch.cat((img, 1 - img))
        img = torch.zeros((1, 64, 64, 64))
        img[0, :10, :10, :10] = 1
        img[0, -10:, -10:, -10:] = 1
        img[0, -10:, :10, -10:] = 1
        data_dict = {'image': torch.cat((img, 1 - img))}
        st = time()
        out = tr(**data_dict)
        times_torch.append(time() - st)
    print('torch', np.mean(times_torch))

    from batchviewer import view_batch
    view_batch(torch.cat((img, 1 - img)), out['image'], torch.cat((img, 1 - img)) - out['image'])