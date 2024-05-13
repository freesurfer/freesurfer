import torch


def crop_tensor(input_tensor, center, crop_size, pad_mode='constant', pad_kwargs=None):
    """
    Crops and pads an input tensor based on the specified center and crop size. Padding can be customized.

    Parameters:
    - input_tensor (torch.Tensor): The input tensor with shape (c, x, y) or (c, x, y, z).
    - center (tuple): The center coordinates of the crop (x, y) or (x, y, z).
    - crop_size (tuple): The size of the crop (width, height) or (width, height, depth).
    - pad_mode (str): The mode to use for padding (see torch.nn.functional.pad documentation).
    - pad_kwargs (dict, optional): Additional keyword arguments for padding.

    Returns:
    - torch.Tensor: The cropped and possibly padded tensor.
    """
    if pad_kwargs is None:
        pad_kwargs = {'value': 0}

    # Calculate dimensions
    dim = len(center)  # Spatial dimensions
    assert len(crop_size) == dim, "Crop size and center must have the same number of dimensions"
    assert input_tensor.ndim - 1 == dim, "Crop size and input_tensor must have the same number of spatial dimensions"

    spatial_shape = input_tensor.shape[-dim:]
    start = [max(0, cen - cs // 2) for cen, cs in zip(center, crop_size)]
    end = [min(sh, st + cs) for sh, st, cs in zip(spatial_shape, start, crop_size)]

    # Calculate padding
    padding_needed = [(cs - (e - s)) for cs, s, e in zip(crop_size, start, end)]
    pad_before = [max(0, - (cen - cs // 2)) for cen, cs in zip(center, crop_size)]
    pad_after = [pn - pb for pn, pb in zip(padding_needed, pad_before)]

    # Adjust start and end for the case where the crop is entirely outside the input tensor
    start = [min(max(0, s), sh) for s, sh in zip(start, spatial_shape)]
    end = [max(min(e, sh), 0) for e, sh in zip(end, spatial_shape)]

    # Perform crop
    slices = [slice(None)] + [slice(s, e) for s, e in zip(start, end)]
    cropped = input_tensor[tuple(slices)]

    # Pad
    pad_width = sum([[b, a] for b, a in zip(pad_before[::-1], pad_after[::-1])], [])
    if any(pad_width):
        cropped = torch.nn.functional.pad(cropped, pad_width, mode=pad_mode, **pad_kwargs)

    return cropped


def center_crop(input_tensor, crop_size, pad_mode='constant', pad_kwargs=None):
    """
    Performs a center crop on the input tensor. If the crop extends beyond the borders of the tensor,
    it will be padded according to the specified pad_mode and pad_kwargs.

    Parameters:
    - input_tensor (torch.Tensor): The input tensor with shape (c, x, y) or (c, x, y, z).
    - crop_size (tuple): The size of the crop (width, height) or (width, height, depth).
    - pad_mode (str): The mode to use for padding (see torch.nn.functional.pad documentation).
    - pad_kwargs (dict, optional): Additional keyword arguments for padding.

    Returns:
    - torch.Tensor: The center-cropped and possibly padded tensor.
    """
    dim = len(input_tensor.shape) - 1  # Number of spatial dimensions (2 or 3)
    spatial_shape = input_tensor.shape[-dim:]  # Spatial dimensions of the input tensor

    # Calculate the center of the input tensor
    center = tuple(s // 2 for s in spatial_shape)

    # Use the previously defined function for cropping and padding
    return crop_tensor(input_tensor, center, crop_size, pad_mode, pad_kwargs)

