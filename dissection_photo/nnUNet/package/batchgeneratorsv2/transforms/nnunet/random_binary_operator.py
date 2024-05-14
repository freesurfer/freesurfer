from time import time
from typing import Union, List, Tuple, Callable

import numpy as np
import torch
from fft_conv_pytorch import fft_conv

from batchgeneratorsv2.helpers.scalar_type import RandomScalar, sample_scalar
from batchgeneratorsv2.transforms.base.basic_transform import ImageOnlyTransform
from skimage.morphology import ball, disk
from skimage.morphology.binary import binary_erosion, binary_dilation, binary_closing, binary_opening

import torch.nn.functional as F


def binary_dilation_torch(input_tensor, structure_element):
    # Convert the boolean tensor to float
    input_tensor = input_tensor.float()

    # Get the number of dimensions of the input tensor
    num_dims = input_tensor.dim()

    # Prepare the structure element for convolution
    # Adding extra dimensions to match the input shape for convolution
    if num_dims == 2:  # For 2D inputs
        structure_element = structure_element.unsqueeze(0).unsqueeze(0).float()
    elif num_dims == 3:  # For 3D inputs, adding batch dimension
        structure_element = structure_element.unsqueeze(0).unsqueeze(0).float()
    else:
        raise ValueError("Input tensor must be 2D (X, Y) or 3D (X, Y, Z).")

    # Perform the convolution
    # if num_dims == 2:  # 2D convolution
    #     output = F.conv2d(input_tensor.unsqueeze(0).unsqueeze(0), structure_element, padding='same')
    # elif num_dims == 3:  # 3D convolution
    #     output = F.conv3d(input_tensor.unsqueeze(0).unsqueeze(0), structure_element, padding='same')
    output = torch.round(fft_conv(input_tensor.unsqueeze(0).unsqueeze(0), structure_element, padding='same'), decimals=0)

    # Threshold to get binary output
    output = output > 0

    # Squeeze the batch dimension out and convert to bool
    return output.squeeze(0).squeeze(0).bool()


def binary_erosion_torch(input_tensor, structure_element):
    return ~binary_dilation_torch(~input_tensor, structure_element)


def binary_opening_torch(input_tensor, structure_element):
    return binary_dilation_torch(binary_erosion_torch(input_tensor, structure_element), structure_element)


def binary_closing_torch(input_tensor, structure_element):
    return binary_erosion_torch(binary_dilation_torch(input_tensor, structure_element), structure_element)


class ApplyRandomBinaryOperatorTransform(ImageOnlyTransform):
    def __init__(self,
                 channel_idx: Union[int, List[int], Tuple[int, ...]],
                 any_of_these: Tuple[Callable, ...] = (binary_dilation_torch, binary_erosion_torch, binary_closing_torch, binary_opening_torch),
                 strel_size: RandomScalar = (1, 10),
                 p_per_label: float = 1):
        """
        We use fft conv. Slower for small kernels but boi does it perform on larger kernels

        Args:
            channel_idx:
            any_of_these:
            strel_size:
            p_per_label:
        """
        super().__init__()
        if not isinstance(channel_idx, (list, tuple)):
            channel_idx = [channel_idx]
        if isinstance(channel_idx, tuple):
            channel_idx = list(channel_idx)
        self.channel_idx = channel_idx
        self.any_of_these = any_of_these
        self.strel_size = strel_size
        self.p_per_label = p_per_label

    def get_parameters(self, **data_dict) -> dict:
        # this needs to be applied in random order to the channels
        np.random.shuffle(self.channel_idx)
        apply_to_channels = [self.channel_idx[i] for i, j in enumerate(torch.rand(len(self.channel_idx)) < self.p_per_label) if j]
        operators = [np.random.choice(self.any_of_these) for _ in apply_to_channels]
        strel_size = [sample_scalar(self.strel_size, image=data_dict['image'], channel=a) for a in apply_to_channels]
        return {
            'apply_to_channels': apply_to_channels,
            'operators': operators,
            'strel_size': strel_size,
        }

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        for a, o, s in zip(params['apply_to_channels'], params['operators'], params['strel_size']):
            # this is a binary map so bool is fine
            workon = img[a]#.numpy()
            orig_dtype = workon.dtype
            workon = workon.to(bool)
            if workon.ndim == 2:
                strel = disk(s, dtype=bool)
            else:
                strel = ball(s, dtype=bool)
            result = o(workon, torch.from_numpy(strel))
            other_ch = [i for i in self.channel_idx if i != a]
            if len(other_ch) > 0:
                was_added_mask = result & (~workon)
                for oc in other_ch:
                    img[oc][was_added_mask] = 0
            img[a] = result.to(orig_dtype)#torch.from_numpy(result)
        return img


# class ApplyRandomBinaryOperatorTransformNpy(ImageOnlyTransform):
#     def __init__(self,
#                  channel_idx: Union[int, List[int], Tuple[int, ...]],
#                  any_of_these: Tuple[Callable, ...] = (binary_dilation, binary_erosion, binary_closing, binary_opening),
#                  strel_size: ScalarType = (1, 10),
#                  p_per_label: float = 1):
#         """
#         We stick to the nnunet implementation and dont optimize. This is a TODO for the future
#
#         Args:
#             channel_idx:
#             any_of_these:
#             strel_size:
#             p_per_label:
#         """
#         super().__init__()
#         if not isinstance(channel_idx, (list, tuple)):
#             channel_idx = [channel_idx]
#         if isinstance(channel_idx, tuple):
#             channel_idx = list(channel_idx)
#         self.channel_idx = channel_idx
#         self.any_of_these = any_of_these
#         self.strel_size = strel_size
#         self.p_per_label = p_per_label
#
#     def get_parameters(self, **data_dict) -> dict:
#         # this needs to be applied in random order to the channels
#         np.random.shuffle(self.channel_idx)
#         apply_to_channels = [self.channel_idx[i] for i, j in enumerate(torch.rand(len(self.channel_idx)) < self.p_per_label) if j]
#         operators = [np.random.choice(self.any_of_these) for _ in apply_to_channels]
#         strel_size = [sample_scalar(self.strel_size, image=data_dict['image'], channel=a) for a in apply_to_channels]
#         return {
#             'apply_to_channels': apply_to_channels,
#             'operators': operators,
#             'strel_size': strel_size,
#         }
#
#     def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
#         for a, o, s in zip(params['apply_to_channels'], params['operators'], params['strel_size']):
#             # this is a binary map so bool is fine
#             workon = img[a].numpy()
#             orig_dtype = workon.dtype
#             workon = workon.astype(bool)
#             if workon.ndim == 2:
#                 strel = disk(s, dtype=bool)
#             else:
#                 strel = ball(s, dtype=bool)
#             result = o(workon, strel)
#             other_ch = [i for i in self.channel_idx if i != a]
#             if len(other_ch) > 0:
#                 was_added_mask = result & (~workon)
#                 for oc in other_ch:
#                     img[oc][was_added_mask] = 0
#             img[a] = torch.from_numpy(result)
#         return img


if __name__ == '__main__':
    torch.set_num_threads(1)

    tr = ApplyRandomBinaryOperatorTransform(
        channel_idx=(0, 1),
        strel_size=1,
    )

    times_torch = []
    for _ in range(10):
        # img = (torch.rand((1, 128, 128, 128)) < 0.5).to(torch.int16)
        # img = torch.cat((img, 1 - img))
        img = torch.zeros((1, 50, 50, 50))
        img[0, :10, :10, :10] = 1
        data_dict = {'image': torch.cat((img, 1 - img))}
        st = time()
        out = tr(**data_dict)
        times_torch.append(time() - st)
    print('torch', np.mean(times_torch))

    # from batchviewer import view_batch
    # view_batch(torch.cat((img, 1 - img)), out['image'], torch.cat((img, 1 - img)) - out['image'])