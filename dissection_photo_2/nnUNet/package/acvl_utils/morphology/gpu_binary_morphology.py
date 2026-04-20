from typing import Union

import numpy as np
import torch
from torch import nn
from dynamic_network_architectures.building_blocks.helper import convert_dim_to_conv_op
from acvl_utils.morphology.morphology_helper import generate_ball
from torch.backends import cudnn


def gpu_binary_dilation(binary_array: Union[np.ndarray, torch.Tensor], selem: np.ndarray) \
        -> Union[np.ndarray, torch.Tensor]:
    """
    IMPORTANT: ALWAYS benchmark your image and kernel sizes first. Sometimes GPU is actually slower than CPU!
    """
    # cudnn.benchmark True DESTROYS our computation time (like 30X decrease lol). Make sure it's disabled (we set is
    # back below)
    benchmark = cudnn.benchmark
    cudnn.benchmark = False
    assert all([i % 2 == 1 for i in selem.shape]), \
        f"Only structure elements of uneven shape supported. Shape is {selem.shape}"

    with torch.no_grad():
        # move source array to GPU first. Uses non-blocking (important!) so that copy operation can run in background.
        # Cast to half only on the GPU because that is much faster and because the source array is quicker to
        # transfger the less bytes per element it has.
        is_tensor = isinstance(binary_array, torch.Tensor)
        if not is_tensor:
            binary_array = torch.from_numpy(binary_array)
        orig_device = binary_array.device
        binary_array = binary_array.to(0, non_blocking=True).half()

        # initialize conv as half
        conv = convert_dim_to_conv_op(len(binary_array.shape))(in_channels=1, out_channels=1,
                                                               kernel_size=selem.shape,
                                                               stride=1,
                                                               padding='same',
                                                               bias=False).half()
        conv.weight = nn.Parameter(torch.from_numpy(selem[None, None]).half())
        conv = conv.to(0, non_blocking=True)

        # no need for autocast because everything is half already (I tried and it doesn't improve computation times)
        # again convert to 1 byte per element byte on GPU, then copy
        out = (conv(binary_array[None, None]) > 0).byte().to(orig_device, non_blocking=True)[0, 0]

    if not is_tensor:
        out = out.numpy()

    torch.cuda.empty_cache()
    # revert changes to cudnn.benchmark
    cudnn.benchmark = benchmark
    return out


def gpu_binary_erosion(binary_array: Union[np.ndarray, torch.Tensor], selem: np.ndarray) \
        -> Union[np.ndarray, torch.Tensor]:
    """
    IMPORTANT: ALWAYS benchmark your image and kernel sizes first. Sometimes GPU is actually slower than CPU!
    """
    # cudnn.benchmark True DESTROYS our computation time (like 30X decrease lol). Make sure it's disabled (we set is
    # back below)
    benchmark = cudnn.benchmark
    cudnn.benchmark = False
    assert all([i % 2 == 1 for i in selem.shape]), \
        f"Only structure elements of uneven shape supported. Shape is {selem.shape}"

    with torch.no_grad():
        # move source array to GPU first. Uses non-blocking (important!) so that copy operation can run in background.
        # Cast to half only on the GPU because that is much faster and because the source array is quicker to
        # transfger the less bytes per element it has.
        is_tensor = isinstance(binary_array, torch.Tensor)
        if not is_tensor:
            binary_array = torch.from_numpy(binary_array)
        orig_device = binary_array.device
        binary_array = binary_array.to(0, non_blocking=True).half()

        # initialize conv as half
        conv = convert_dim_to_conv_op(len(binary_array.shape))(in_channels=1, out_channels=1,
                                                               kernel_size=selem.shape,
                                                               stride=1,
                                                               padding='same',
                                                               bias=False).half()
        conv.weight = nn.Parameter(torch.from_numpy(selem[None, None]).half())
        conv = conv.to(0, non_blocking=True)

        # no need for autocast because everything is half already (I tried and it doesn't improve computation times)
        # again convert to 1 byte per element byte on GPU, then copy
        out = (conv(binary_array[None, None]) == selem.sum()).to(orig_device)[0, 0]

    if not is_tensor:
        out = out.numpy()
    torch.cuda.empty_cache()
    # revert changes to cudnn.benchmark
    cudnn.benchmark = benchmark
    return out


def gpu_binary_opening(binary_array: Union[np.ndarray, torch.Tensor], selem: np.ndarray) -> Union[np.ndarray, torch.Tensor]:
    return gpu_binary_dilation(gpu_binary_erosion(binary_array, selem), selem)


def gpu_binary_closing(binary_array: Union[np.ndarray, torch.Tensor], selem: np.ndarray) -> Union[np.ndarray, torch.Tensor]:
    return gpu_binary_erosion(gpu_binary_dilation(binary_array, selem), selem)


if __name__ == '__main__':
    from skimage.morphology import binary_dilation, binary_erosion
    cudnn.benchmark = True
    from time import time
    # inp = np.zeros((512, 512, 512), dtype=bool)
    # inp[10:100, 50:100, 100:120] = True
    # selem = generate_ball((7, 7, 7))

    inp = np.zeros((160, 2560, 2160), dtype=bool)
    inp[200:300, 200:300, 100:200] = True
    selem = generate_ball((5, 5, 3))

    # one dry run for warmup
    _ = gpu_binary_dilation(inp, selem)

    start = time()
    output_gpu = gpu_binary_dilation(inp, selem)
    time_gpu = time() - start
    print(f'Dilation: GPU: {time_gpu}s')

    start = time()
    ref = binary_dilation(inp, selem)
    time_skimage = time() - start
    print(f'Dilation: CPU: {time_skimage}s')

    assert np.all(output_gpu == ref)

    start = time()
    output_gpu = gpu_binary_erosion(inp, selem)
    time_gpu = time() - start
    print(f'Erosion: GPU: {time_gpu}s')

    start = time()
    ref = binary_erosion(inp, selem)
    time_skimage = time() - start
    print(f'Erosion: CPU: {time_skimage}s')

    assert np.all(output_gpu == ref)


