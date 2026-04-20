from typing import Tuple

import torch
from torch.nn import functional as F


def maybe_resample_on_gpu(data: torch.Tensor,
                          target_shape: Tuple[int, ...],
                          return_type=torch.half,
                          compute_device: str = 'cuda:0', result_device: str = 'cuda:0',
                          compute_precision=torch.half,
                          fallback_compute_device: str = 'cpu', fallback_result_device: str = 'cpu',
                          fallback_compute_precision=torch.float,
                          verbose: bool = True) -> torch.Tensor:
    """
    channel-wise resampling function. Most useful for softmax resampling. Tries three things:
    - resample the entire array at once on compute_device (fastest but uses the most memory). returned tensor will be on result_device
    - failing that it will resample channel-wise on compute_device. returned tensor will be on fallback_result_device
    - failing that it will resample channel-wise on fallback_device. returned tensor will be on fallback_result_device

    data must be a 4d tensor of shape (c, x, y, z). It will be resampled to shape (tx, ty, tz), so the returned tensor
    will have shape (c, tx, ty, tz)

    returned tensor will be on result_device/fallback_result_device
    returned tensor will have dtype return_type. half recommended if memory is scarce

    fallback_compute_precision must be torch.float for CPU (half not implemented as of 11/2022)

    this function does not use with torch.no_grad()!! You need to do this externally

    Defaults for compute_device, result_device, compute_precision, fallback_compute_device, fallback_result_device,
    fallback_compute_precision are good and should only be deviated from by experts.
    """

    def _maybe_empty_cache():
        if torch.cuda.is_available(): torch.cuda.empty_cache()

    def _allocate_result_tensor(shape, dtype, device, fallback_device) -> torch.Tensor:
        try:
            softmax_resampled = torch.zeros(shape, dtype=dtype, device=device)
        except RuntimeError:
            if verbose:
                print(f'could not allocate result tensor of shape {(data.shape[0], *target_shape)} on device '
                      f'{result_device}. Using fallback {fallback_device}')
            assert fallback_device != device, 'Allocation on device failed and fallback_device=device so there is no ' \
                                              'point attempting this. Duh.'
            softmax_resampled = torch.zeros(shape, dtype=dtype, device=fallback_device)
        finally:
            _maybe_empty_cache()
        return softmax_resampled

    _maybe_empty_cache()

    if all([i == j for i, j in zip(data.shape[1:], target_shape)]):
        try:
            return data.to(result_device)
        except RuntimeError:
            return data.to(fallback_result_device)

    try:
        # try all in one
        result_tensor = \
        F.interpolate(data.to(compute_device).type(compute_precision)[None], target_shape, mode='trilinear')[0].type(
            return_type).to(result_device)
    except RuntimeError:
        if verbose:
            print(f'Resampling all in one failed, attempting channel-wise on compute_device {compute_device}')
        try:
            _maybe_empty_cache()
            # allocate result tensor
            result_tensor = _allocate_result_tensor((data.shape[0], *target_shape), return_type, result_device,
                                                    fallback_result_device)
            for c in range(data.shape[0]):
                result_tensor[c] = F.interpolate(data[c][None, None].to(compute_device).type(compute_precision),
                                                 size=target_shape,
                                                 mode='trilinear')[0, 0].to(result_tensor.device)
        except RuntimeError:
            if verbose:
                print(f'Resampling channel-wise on compute_device {compute_device} failed, attempting chanel-wise '
                      f'on fallback {fallback_compute_device}')
            _maybe_empty_cache()
            # allocate result tensor
            result_tensor = _allocate_result_tensor((data.shape[0], *target_shape), return_type, fallback_result_device,
                                                    fallback_result_device)
            for c in range(data.shape[0]):
                result_tensor[c] = F.interpolate(data[c][None, None].to(fallback_compute_device).type(fallback_compute_precision),
                                                 size=target_shape,
                                                 mode='trilinear')[0, 0].to(result_tensor.device)
        finally:
            _maybe_empty_cache()
    finally:
        _maybe_empty_cache()
    return result_tensor
