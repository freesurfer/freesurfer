import torch
from tqdm import tqdm
import SimpleITK as sitk
import numpy as np
from typing import Any, Union, List, TypeVar, Optional
from time import time
import torch
from torch import nn
from torch.optim import SGD, Adam
from torch.nn import functional as F
from functools import partial
from copy import deepcopy
from abc import ABC, abstractmethod
from torch.optim import SGD, Adam
from collections import deque
from torch.utils.checkpoint import checkpoint


T = TypeVar("T")
devicetype = Union[str, torch.device]
ItemOrList = Union[T, List[T]]


MIN_IMG_SIZE = 32

@torch.jit.script
def make_rectangular_kernel(kernel_size: int) -> torch.Tensor:
    return torch.ones(kernel_size)

@torch.jit.script
def make_triangular_kernel(kernel_size: int) -> torch.Tensor:
    fsize = (kernel_size + 1) // 2
    if fsize % 2 == 0:
        fsize -= 1
    f = torch.ones((1, 1, fsize), dtype=torch.float).div(fsize)
    padding = (kernel_size - fsize) // 2 + fsize // 2
    return F.conv1d(f, f, padding=padding).reshape(-1)

@torch.jit.script
def gaussian_1d(
    sigma: torch.Tensor, truncated: float = 4.0, approx: str = "erf", normalize: bool = True
) -> torch.Tensor:
    sigma = torch.as_tensor(sigma, dtype=torch.float, device=sigma.device if isinstance(sigma, torch.Tensor) else None)
    device = sigma.device
    if truncated <= 0.0:
        raise ValueError(f"truncated must be positive, got {truncated}.")
    tail = int(max(float(sigma) * truncated, 0.5) + 0.5)
    if approx.lower() == "erf":
        x = torch.arange(-tail, tail + 1, dtype=torch.float, device=device)
        t = 0.70710678 / torch.abs(sigma)
        out = 0.5 * ((t * (x + 0.5)).erf() - (t * (x - 0.5)).erf())
        out = out.clamp(min=0)
    elif approx.lower() == "sampled":
        x = torch.arange(-tail, tail + 1, dtype=torch.float, device=sigma.device)
        out = torch.exp(-0.5 / (sigma * sigma) * x**2)
        if not normalize:  # compute the normalizer
            out = out / (2.5066282 * sigma)
    else:
        raise NotImplementedError(f"Unsupported option: approx='{approx}'.")
    return out / out.sum() if normalize else out  # type: ignore


@torch.jit.script
def make_gaussian_kernel(kernel_size: int) -> torch.Tensor:
    sigma = torch.tensor(kernel_size / 3.0)
    kernel = gaussian_1d(sigma=sigma, truncated=(kernel_size // 2) * 1.0, approx="sampled", normalize=False) * (
        2.5066282 * sigma
    )
    return kernel[:kernel_size]

kernel_dict = {
    "rectangular": make_rectangular_kernel,
    "triangular": make_triangular_kernel,
    "gaussian": make_gaussian_kernel,
}


def integer_to_onehot(image: torch.Tensor, background_label:int=0, max_label=None):
    ''' convert an integer map into one hot mapping
    assumed the image is of size [H, W, [D]] and we convert it into [C, H, W, [D]]

    background_label: this is the label to ignore (default: 0)
    max_label: max value of the label expected in the label segmentation, which sometimes may not be present
    we provide this as an additional option in case some images do not have the anatomy corresponding to the max label

    if None, we assume the image has that label already
    '''
    if max_label is None:
        max_label = image.max()
    if background_label >= 0 and background_label <= max_label: # we will ignore it
        num_labels = max_label
    else:
        num_labels = max_label + 1
    onehot = torch.zeros((num_labels, *image.shape), dtype=torch.float32, device=image.device)
    count = 0
    for i in range(num_labels+1):
        if i == background_label:
            continue
        onehot[count, ...] = (image == i)
        count += 1
    return onehot

def _assert_check_scales_decreasing(scales: List[int]):
    ''' Check if the list of scales is in decreasing order '''
    for i in range(len(scales)-1):
        if scales[i] <= scales[i+1]:
            raise ValueError("Scales must be in decreasing order")

def jacobian_2d(u: torch.Tensor, normalize: bool):
    ''' u: displacement vector of size [N, H, W, 2] '''
    B, H, W, _ = u.shape
    newshape = [B, 2, H, W, 2]
    J = torch.empty(newshape, dtype=u.dtype, device=u.device)
    # Compute Jacobian of u and v using image_gradient_singlechannel function
    for i in range(2):
        J[..., i] = image_gradient_singlechannel(u[..., i].reshape(B, 1, H, W), normalize)
    return J

def jacobian_3d(u: torch.Tensor, normalize: bool):
    ''' u: displacement vector of size [N, H, W, D, 3] '''
    B, H, W, D, _ = u.shape
    newshape = [B, 3, H, W, D, 3]
    J = torch.empty(newshape, dtype=u.dtype, device=u.device)
    for i in range(3):
        J[..., i] = image_gradient_singlechannel(u[..., i].reshape(B, 1, H, W, D), normalize)
    return J

def jacobian_fn(u: torch.Tensor, normalize=True):
    '''
    u: displacement vector of size [N, H, W, [D], dims]
    '''
    if len(u.shape) == 4:
        return jacobian_2d(u, normalize)
    elif len(u.shape) == 5:
        return jacobian_3d(u, normalize)
    else:
        raise ValueError(f"jacobian not implemented for tensors of shape {u.shape}")


@torch.jit.script
def _separable_filtering_conv(
    input_: torch.Tensor,
    kernels: List[torch.Tensor],
    pad_mode: str,
    spatial_dims: int,
    paddings: List[int],
    num_channels: int,
) -> torch.Tensor:

    # re-write from recursive to non-recursive for torch.jit to work
    # for d in range(spatial_dims-1, -1, -1):
    for d in range(spatial_dims):
        s = [1] * len(input_.shape)
        s[d + 2] = -1
        _kernel = kernels[d].reshape(s)
        # if filter kernel is unity, don't convolve
        if _kernel.numel() == 1 and _kernel[0] == 1:
            continue

        _kernel = _kernel.repeat([num_channels, 1] + [1] * spatial_dims)
        _padding = [0] * spatial_dims
        _padding[d] = paddings[d]
        _reversed_padding = _padding[::-1]

        # translate padding for input to torch.nn.functional.pad
        _reversed_padding_repeated_twice: list[list[int]] = [[p, p] for p in _reversed_padding]
        _sum_reversed_padding_repeated_twice: list[int] = []
        for p in _reversed_padding_repeated_twice:
            _sum_reversed_padding_repeated_twice.extend(p)
        # _sum_reversed_padding_repeated_twice: list[int] = sum(_reversed_padding_repeated_twice, [])

        padded_input = F.pad(input_, _sum_reversed_padding_repeated_twice, mode=pad_mode)
        # update input
        if spatial_dims == 1:
            input_ = F.conv1d(input=padded_input, weight=_kernel, groups=num_channels)
        elif spatial_dims == 2:
            input_ = F.conv2d(input=padded_input, weight=_kernel, groups=num_channels)
        elif spatial_dims == 3:
            input_ = F.conv3d(input=padded_input, weight=_kernel, groups=num_channels)
        else:
            raise NotImplementedError(f"Unsupported spatial_dims: {spatial_dims}.")
    return input_

@torch.jit.script
def separable_filtering(x: torch.Tensor, kernels: ItemOrList[torch.Tensor], mode: str = "zeros") -> torch.Tensor:
    if not isinstance(x, torch.Tensor):
        raise TypeError(f"x must be a torch.Tensor but is {type(x).__name__}.")

    spatial_dims = len(x.shape) - 2
    if isinstance(kernels, torch.Tensor):
        kernels = [kernels] * spatial_dims
    _kernels = [s.to(x) for s in kernels]
    _paddings = [(k.shape[0] - 1) // 2 for k in _kernels]
    n_chs = x.shape[1]
    pad_mode = "constant" if mode == "zeros" else mode
    return _separable_filtering_conv(x, _kernels, pad_mode, spatial_dims, _paddings, n_chs)

def scaling_and_squaring(u, grid, n=6):
    """
    Apply scaling and squaring to a displacement field

    :param u: Input stationary velocity field, PyTorch tensor of shape [B, D, H, W, 3] or [B, H, W, 2]
    :param grid: Sampling grid of size [B, D, H, W, dims]  or [B, H, W, dims]
    :param n: Number of iterations of scaling and squaring (default: 6)

    :returns: Output displacement field, v, PyTorch tensor of shape [B, D, H, W, dims] or [B, H, W, dims]
    """
    dims = u.shape[-1]
    v = (1.0 / 2 ** n) * u
    if dims == 3:
        for i in range(n):
            vimg = v.permute(0, 4, 1, 2, 3)  # [1, 3, D, H, W]
            v = v + F.grid_sample(vimg, v + grid, align_corners=True).permute(0, 2, 3, 4, 1)
    elif dims == 2:
        for i in range(n):
            vimg = v.permute(0, 3, 1, 2)
            v = v + F.grid_sample(vimg, v + grid, align_corners=True).permute(0, 2, 3, 1)
    else:
        raise ValueError('Invalid dimension: {}'.format(dims))
    return v

def compute_inverse_warp_exp(warp, grid, lr=5e-3, iters=200, n=10):
    ''' compute warp inverse using exponential map '''
    with torch.set_grad_enabled(True):
        vel = nn.Parameter(torch.zeros_like(warp))
        optim = torch.optim.Adam([vel], lr=lr)
        permute_vtoimg = (0, 4, 1, 2, 3) if len(warp.shape) == 5 else (0, 3, 1, 2)
        permute_imgtov = (0, 2, 3, 4, 1) if len(warp.shape) == 5 else (0, 2, 3, 1)
        # pbar = tqdm(range(iters))
        pbar = range(iters)
        for i in pbar:
            optim.zero_grad()
            invwarp = scaling_and_squaring(vel, grid, n=n)
            loss = invwarp + F.grid_sample(warp.permute(*permute_vtoimg), grid + invwarp, mode='bilinear', align_corners=True).permute(*permute_imgtov)
            loss2 = warp + F.grid_sample(invwarp.permute(*permute_vtoimg), grid + warp, mode='bilinear', align_corners=True).permute(*permute_imgtov)
            loss = (loss**2).sum() + (loss2**2).sum()
            loss.backward()
            optim.step()
    return scaling_and_squaring(vel.data, grid, n=n)

def grad_smoothing_hook(grad: torch.Tensor, gaussians: List[torch.Tensor]):
    ''' this backward hook will smooth out the gradient using the gaussians
    has to be called with a partial function
    '''
    # grad is of shape [B, H, W, D, dims]
    if len(grad.shape) == 5:
        permute_vtoimg = (0, 4, 1, 2, 3)
        permute_imgtov = (0, 2, 3, 4, 1)
    elif len(grad.shape) == 4:
        permute_vtoimg = (0, 3, 1, 2)
        permute_imgtov = (0, 2, 3, 1)
    return separable_filtering(grad.permute(*permute_vtoimg), gaussians).permute(*permute_imgtov)

def downsample(image: ItemOrList[torch.Tensor], size: List[int], mode: str, sigma: Optional[torch.Tensor]=None,
               gaussians: Optional[torch.Tensor] = None) -> torch.Tensor:
    '''
    this function is to downsample the image to the given size
    but first, we need to perform smoothing
    if sigma is provided (in voxels), then use this sigma for downsampling, otherwise infer sigma
    '''
    if gaussians is None:
        if sigma is None:
            orig_size = list(image.shape[2:])
            sigma = [0.5 * orig_size[i] / size[i] for i in range(len(orig_size))]   # use sigma as the downsampling factor
        sigma = torch.tensor(sigma, dtype=torch.float32, device=image.device)
        # create gaussian convs
        gaussians = [gaussian_1d(s, truncated=2) for s in sigma]
    # otherwise gaussians is given, just downsample
    image_smooth = separable_filtering(image, gaussians)
    image_down = F.interpolate(image_smooth, size=size, mode=mode, align_corners=True)
    return image_down

class Image:
    '''
    TODO: Documentation here
    '''

    def __init__(self, itk_image: sitk.SimpleITK.Image, device: devicetype = 'cuda',
                 is_segmentation=False, max_seg_label=None, background_seg_label=0, seg_preprocessor=lambda x: x,
                 spacing=None, direction=None, origin=None, center=None) -> None:
        self.itk_image = itk_image
        # check for segmentation parameters
        # if `is_segmentation` is False, then just treat this as a float image
        if not is_segmentation:
            self.array = torch.from_numpy(sitk.GetArrayFromImage(itk_image).astype(float)).to(device).float()
            self.array = self.array[
                None, None]  # TODO: Change it to support multichannel images, right now just batchify and add a dummy channel to it
            channels = itk_image.GetNumberOfComponentsPerPixel()
            self.channels = channels
            assert channels == 1, "Only single channel images supported"
        else:
            array = torch.from_numpy(sitk.GetArrayFromImage(itk_image).astype(int)).to(device).long()
            # preprocess segmentation if provided by user
            array = seg_preprocessor(array)
            if max_seg_label is not None:
                array[array > max_seg_label] = background_seg_label
            array = integer_to_onehot(array, background_label=background_seg_label, max_label=max_seg_label)[None]  # []
            self.array = array.float()
            self.channels = array.shape[1]
        # initialize matrix for pixel to physical
        dims = itk_image.GetDimension()
        self.dims = dims
        if dims not in [2, 3]:
            raise NotImplementedError("Image class only supports 2D/3D images.")

        # custom spacing if not provided use simpleitk values
        spacing = np.array(itk_image.GetSpacing())[None] if spacing is None else np.array(spacing)[None]
        origin = np.array(itk_image.GetOrigin())[None] if origin is None else np.array(origin)[None]
        direction = np.array(itk_image.GetDirection()).reshape(dims, dims) if direction is None else np.array(
            direction).reshape(dims, dims)
        if center is not None:
            print("Center location provided, recalibrating origin.")
            origin = center - np.matmul(direction, ((np.array(itk_image.GetSize()) * spacing / 2).squeeze())[:, None]).T

        px2phy = np.eye(dims + 1)
        px2phy[:dims, -1] = origin
        px2phy[:dims, :dims] = direction
        px2phy[:dims, :dims] = px2phy[:dims, :dims] * spacing
        # generate mapping from torch to px
        torch2px = np.eye(dims + 1)
        scaleterm = (np.array(itk_image.GetSize()) - 1) * 0.5
        torch2px[:dims, :dims] = np.diag(scaleterm)
        torch2px[:dims, -1] = scaleterm
        # save the mapping from physical to torch and vice versa
        self.torch2phy = torch.from_numpy(np.matmul(px2phy, torch2px)).to(device).float().unsqueeze(0)
        self.phy2torch = torch.inverse(self.torch2phy[0]).float().unsqueeze(0)
        # also save intermediates just in case (as numpy arrays)
        self._torch2px = torch2px
        self._px2phy = px2phy
        self.device = device

    @classmethod
    def load_file(cls, image_path: str, *args, **kwargs) -> 'Image':
        itk_image = sitk.ReadImage(image_path)
        return cls(itk_image, *args, **kwargs)

    @property
    def shape(self):
        return self.array.shape


class BatchedImages:
    '''
    Class for batched images
    '''

    def __init__(self, images: Union[Image, List[Image]]) -> None:
        if isinstance(images, Image):
            images = [images]
        self.images = images
        if len(self.images) == 0:
            raise ValueError("BatchedImages must have at least one image")
        for image in self.images:
            if not isinstance(image, Image):
                raise TypeError("All images must be of type Image")
        shapes = [x.array.shape for x in self.images]
        if all([x == shapes[0] for x in shapes]):
            self.shape = shapes[0]
        else:
            raise ValueError("All images must have the same shape")
        self.n_images = len(self.images)
        self.interpolate_mode = 'bilinear' if len(self.images[0].shape) == 4 else 'trilinear'

    def __call__(self):
        # get batch of images
        return torch.cat([x.array for x in self.images], dim=0)

    @property
    def device(self):
        return self.images[0].device

    @property
    def dims(self):
        return self.images[0].dims

    def size(self):
        return self.n_images

    def shape(self):
        shape = self.images[0].shape
        shape[0] = self.n_images
        return shape

    def get_torch2phy(self):
        return torch.cat([x.torch2phy for x in self.images], dim=0)

    def get_phy2torch(self):
        return torch.cat([x.phy2torch for x in self.images], dim=0)



class AbstractDeformation(ABC):
    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def get_warp(self):
        ''' returns displacement field '''
        raise NotImplementedError

    @abstractmethod
    def get_inverse_warp(self):
        ''' returns inverse displacement field '''
        raise NotImplementedError

    @abstractmethod
    def set_size(self, size):
        ''' sets size of the deformation field '''
        raise NotImplementedError

    @abstractmethod
    def step(self):
        ''' optimizes the deformation field '''
        raise NotImplementedError


class ConvergenceMonitor:
    def __init__(self, N, slope):
        """
        Initialize the ConvergenceMonitor class.
        Args:
        - N: number of values to keep track of.
        """
        self.N = N
        self.losses = deque(maxlen=N)
        self.slope = slope

    def update(self, loss):
        """Append a new loss value to the monitor."""
        self.losses.append(loss)

    def _compute_slope(self):
        """Compute the slope of the best-fit line using simple linear regression."""
        if len(self.losses) < 2:
            # Can't compute a slope with less than 2 points
            return 0

        x = np.arange(len(self.losses))
        y = np.array(self.losses)

        # Compute the slope (m) of the best-fit line y = mx + c
        # m = (NΣxy - ΣxΣy) / (NΣx^2 - (Σx)^2)
        xy_sum = np.dot(x, y)
        x_sum = x.sum()
        y_sum = y.sum()
        x_squared_sum = (x ** 2).sum()
        N = len(self.losses)

        numerator = N * xy_sum - x_sum * y_sum
        denominator = N * x_squared_sum - x_sum ** 2
        if denominator == 0:
            return 0

        slope = numerator / denominator
        return slope

    def converged(self, loss=None):
        """Check if the loss has increased (i.e., slope > threshold).
        optionally, update the monitor with a new loss value.
        """
        if loss is not None:
            self.update(loss)
        if len(self.losses) < self.N:
            return False
        else:
            slope = self._compute_slope()
            return slope > self.slope

    def reset(self):
        self.losses.clear()


class AbstractRegistration(ABC):

    def __init__(self,
                 scales: List[int], iterations: List[float],
                 fixed_images: BatchedImages, moving_images: BatchedImages,
                 loss_type: str = "cc",
                 mi_kernel_type: str = 'b-spline', cc_kernel_type: str = 'rectangular',
                 custom_loss: nn.Module = None,
                 loss_params: dict = {},
                 cc_kernel_size: int = 3,
                 reduction: str = 'mean',
                 tolerance: float = 1e-6, max_tolerance_iters: int = 10,
                 progress_bar: bool = True,
                 ) -> None:
        '''
        Initialize abstract registration class
        '''
        super().__init__()
        self.scales = scales
        _assert_check_scales_decreasing(self.scales)
        self.iterations = iterations
        assert len(self.iterations) == len(self.scales), "Number of iterations must match number of scales"
        # check for fixed and moving image sizes
        self.fixed_images = fixed_images
        self.moving_images = moving_images
        assert (self.fixed_images.size() == self.moving_images.size()), "Number of fixed and moving images must match"

        self.tolerance = tolerance
        self.max_tolerance_iters = max_tolerance_iters
        self.convergence_monitor = ConvergenceMonitor(self.max_tolerance_iters, self.tolerance)

        self.device = fixed_images.device
        self.dims = self.fixed_images.dims
        self.progress_bar = progress_bar  # variable to show or hide progress bar

        # initialize losses
        if loss_type == 'mi':
            raise ValueError('Eugenio lost this LNCC loss when trimming the package')
            # self.loss_fn = GlobalMutualInformationLoss(kernel_type=mi_kernel_type, reduction=reduction, **loss_params)
        elif loss_type == 'cc':
            raise ValueError('Eugenio lost this LNCC loss when trimming the package')
            # self.loss_fn = LocalNormalizedCrossCorrelationLoss(kernel_type=cc_kernel_type, spatial_dims=self.dims,
            #                                                    kernel_size=cc_kernel_size, reduction=reduction,
            #                                                    **loss_params)
        elif loss_type == 'custom':
            self.loss_fn = custom_loss
        elif loss_type == 'mse':
            raise Exception('Eugenio lost this loss when trimming the package')
            # self.loss_fn = partial(F.mse_loss, reduction=reduction)
        else:
            raise ValueError(f"Loss type {loss_type} not supported")

    @abstractmethod
    def optimize(self):
        pass

class GeodesicShooting(nn.Module, AbstractDeformation):
    '''
    Class for geodesic shooting, by optimizing a velocity field
    '''

    def __init__(self, fixed_images: BatchedImages, moving_images: BatchedImages,
                 integrator_n: Union[str, int] = 6,
                 optimizer: str = 'Adam', optimizer_lr: float = 1e-2, optimizer_params: dict = {},
                 smoothing_grad_sigma: float = 0.5,
                 init_scale: int = 1,
                 ) -> None:
        super().__init__()
        self.num_images = num_images = fixed_images.size()
        spatial_dims = fixed_images.shape[2:]  # [H, W, [D]]
        if init_scale > 1:
            spatial_dims = [max(int(s / init_scale), MIN_IMG_SIZE) for s in spatial_dims]
        self.n_dims = len(spatial_dims)  # number of spatial dimensions
        self.device = fixed_images.device
        # permute indices  (image to v and v to image)
        self.permute_imgtov = (0, *range(2, self.n_dims + 2), 1)  # [N, HWD, dims] -> [N, HWD, dims] -> [N, dims, HWD]
        self.permute_vtoimg = (0, self.n_dims + 1, *range(1, self.n_dims + 1))  # [N, dims, HWD] -> [N, HWD, dims]
        # define velocity field
        velocity_field = torch.zeros([num_images, *spatial_dims, self.n_dims], dtype=torch.float32,
                                     device=fixed_images.device)  # [N, HWD, dims]
        # attach grad hook if smoothing is required
        self.smoothing_grad_sigma = smoothing_grad_sigma
        if smoothing_grad_sigma > 0:
            self.smoothing_grad_gaussians = [gaussian_1d(s, truncated=2) for s in (
                        torch.zeros(self.n_dims, device=fixed_images.device) + smoothing_grad_sigma)]
        # init grid, velocity field and grad hook
        self.initialize_grid(spatial_dims)
        self.register_parameter('velocity_field', nn.Parameter(velocity_field))
        self.attach_grad_hook()
        # self.velocity_field = nn.Parameter(velocity_field)
        self.integrator_n = integrator_n
        # define optimizer
        self.optimizer = getattr(torch.optim, optimizer)([self.velocity_field], lr=optimizer_lr,
                                                         **deepcopy(optimizer_params))
        self.optimizer_lr = optimizer_lr
        self.optimizer_name = optimizer

    def attach_grad_hook(self):
        ''' attack the grad hook to the velocity field if needed '''
        if self.smoothing_grad_sigma > 0:
            hook = partial(grad_smoothing_hook, gaussians=self.smoothing_grad_gaussians)
            self.velocity_field.register_hook(hook)

    def initialize_grid(self, size):
        ''' initialize grid to a size '''
        grid = F.affine_grid(
            torch.eye(self.n_dims, self.n_dims + 1, device=self.device)[None].expand(self.num_images, -1, -1), \
            [self.num_images, self.n_dims, *size], align_corners=True)
        self.register_buffer('grid', grid)
        # self.grid = grid

    def set_zero_grad(self):
        self.optimizer.zero_grad()

    def step(self):
        self.optimizer.step()

    def get_warp(self):
        ''' integrate the velocity field to get the warp '''
        if self.integrator_n == 'auto':
            raise Exception('not implemented (Eugenio: it truly isnt!)')
            # n = _find_integrator_n(self.velocity_field)
        else:
            n = self.integrator_n
        warp = scaling_and_squaring(self.velocity_field, self.grid, n=n)
        return warp

    def get_inverse_warp(self, *args, **kwargs):
        # if self.integrator_n == 'auto':
        #     n = _find_integrator_n(self.velocity_field)
        # else:
        #     n = self.integrator_n
        # invwarp = scaling_and_squaring(-self.velocity_field, self.grid, n=n)
        # return invwarp
        return compute_inverse_warp_exp(self.get_warp().detach(), self.grid)

    def set_size(self, size):
        ''' size: [H, W, D] or [H, W] '''
        mode = 'bilinear' if self.n_dims == 2 else 'trilinear'
        # keep old items for copying
        old_shape = self.velocity_field.shape
        old_optimizer_state = self.optimizer.state_dict()
        # get new velocity field
        velocity_field = F.interpolate(self.velocity_field.detach().permute(*self.permute_vtoimg), size=size, mode=mode,
                                       align_corners=True,
                                       ).permute(*self.permute_imgtov)
        velocity_field = nn.Parameter(velocity_field)
        self.register_parameter('velocity_field', velocity_field)
        self.attach_grad_hook()
        # self.velocity_field = velocity_field
        self.initialize_grid(size)
        self.optimizer = getattr(torch.optim, self.optimizer_name)([self.velocity_field], lr=self.optimizer_lr)
        # TODO: copy state variables from old optimizer
        state_dict = old_optimizer_state['state']
        old_optimizer_state['param_groups'] = self.optimizer.state_dict()['param_groups']
        for g in state_dict.keys():
            for k, v in state_dict[g].items():
                # this is probably a state of the tensor
                if isinstance(v, torch.Tensor) and v.shape == old_shape:
                    state_dict[g][k] = F.interpolate(v.permute(*self.permute_vtoimg), size=size, mode=mode,
                                                     align_corners=True,
                                                     ).permute(*self.permute_imgtov)
        #         if isinstance(v, torch.Tensor):
        #             print(k, v.shape)
        #         else:
        #             print(k, v)
        # input("Here.")
        self.optimizer.load_state_dict(old_optimizer_state)


class WarpAdam:
    ''' at the moment we only support a single warp function
    also supports multi-scale (by simply interpolating to the target size)
    shape of warp = [B, H, W, [D], dims]
    '''

    def __init__(self, warp, lr, warpinv=None, beta1=0.9, beta2=0.99, weight_decay=0, eps=1e-8,
                 scaledown=False, multiply_jacobian=False,
                 smoothing_gaussians=None, optimize_inverse_warp=False):
        # init
        if beta1 < 0.0 or beta1 >= 1.0:
            raise ValueError("Invalid beta1 value: {}".format(beta1))
        if beta2 < 0.0 or beta2 >= 1.0:
            raise ValueError("Invalid beta2 value: {}".format(beta2))
        if weight_decay < 0.0:
            raise ValueError("Invalid weight_decay value: {}".format(weight_decay))
        if lr < 0.0:
            raise ValueError("Invalid lr value: {}".format(lr))
        self.n_dims = len(warp.shape) - 2
        # get half resolutions
        self.half_resolution = 1.0 / (max(warp.shape[1:-1]) - 1)
        self.warp = warp
        self.warpinv = warpinv
        self.optimize_inverse_warp = optimize_inverse_warp
        self.lr = lr
        self.eps = eps
        self.beta1 = beta1
        self.beta2 = beta2
        self.step_t = 0  # initialize step to 0
        self.weight_decay = weight_decay
        self.multiply_jacobian = multiply_jacobian
        self.scaledown = scaledown  # if true, the scale the gradient even if norm is below 1
        self.exp_avg = torch.zeros_like(warp)
        self.exp_avg_sq = torch.zeros_like(warp)
        self.permute_imgtov = (0, *range(2, self.n_dims + 2), 1)  # [N, HWD, dims] -> [N, HWD, dims] -> [N, dims, HWD]
        self.permute_vtoimg = (0, self.n_dims + 1, *range(1, self.n_dims + 1))  # [N, dims, HWD] -> [N, HWD, dims]
        # set grid
        self.batch_size = batch_size = warp.shape[0]
        # init grid
        self.affine_init = torch.eye(self.n_dims, self.n_dims + 1, device=warp.device)[None].expand(batch_size, -1, -1)
        self.initialize_grid(warp.shape[1:-1])
        # gaussian smoothing parameters (if any)
        self.smoothing_gaussians = smoothing_gaussians

    def set_data_and_size(self, warp, size, grid_copy=None, warpinv=None):
        ''' change the optimization variables sizes '''
        self.warp = warp
        mode = 'bilinear' if self.n_dims == 2 else 'trilinear'
        self.exp_avg = F.interpolate(self.exp_avg.detach().permute(*self.permute_vtoimg), size=size, mode=mode,
                                     align_corners=True,
                                     ).permute(*self.permute_imgtov)
        self.exp_avg_sq = F.interpolate(self.exp_avg_sq.detach().permute(*self.permute_vtoimg), size=size, mode=mode,
                                        align_corners=True,
                                        ).permute(*self.permute_imgtov)
        self.half_resolution = 1.0 / (max(warp.shape[1:-1]) - 1)
        self.initialize_grid(size, grid_copy=grid_copy)
        # print(self.warp.shape, warpinv)
        if self.optimize_inverse_warp and warpinv is not None:
            self.warpinv = warpinv

    def initialize_grid(self, size, grid_copy=None):
        ''' initialize the grid (so that we can use it independent of the grid elsewhere) '''
        if grid_copy is None:
            self.grid = F.affine_grid(self.affine_init, [self.batch_size, 1, *size], align_corners=True).detach()
        else:
            self.grid = grid_copy

    def zero_grad(self):
        ''' set the gradient to none '''
        self.warp.grad = None

    def augment_jacobian(self, u):
        # Multiply u (which represents dL/dphi most likely) with Jacobian indexed by J[..., xyz, ..., phi]
        jac = jacobian_fn(self.warp.data + self.grid, normalize=True)  # [B, dims, H, W, [D], dims]
        if self.n_dims == 2:
            ujac = torch.einsum('bxhwp,bhwp->bhwx', jac, u)
        else:
            ujac = torch.einsum('bxhwdp,bhwdp->bhwdx', jac, u)
        return ujac

    def step(self):
        ''' check for momentum, and other things '''
        grad = torch.clone(self.warp.grad.data).detach()
        if self.multiply_jacobian:
            grad = self.augment_jacobian(grad)
        # add weight decay term
        if self.weight_decay > 0:
            grad.add_(self.warp.data, alpha=self.weight_decay)
        # compute moments
        self.step_t += 1
        self.exp_avg.mul_(self.beta1).add_(grad, alpha=1 - self.beta1)
        self.exp_avg_sq.mul_(self.beta2).addcmul_(grad, grad.conj(), value=1 - self.beta2)
        # bias correction
        beta_correction1 = 1 - self.beta1 ** self.step_t
        beta_correction2 = 1 - self.beta2 ** self.step_t
        denom = (self.exp_avg_sq / beta_correction2).sqrt().add_(self.eps)
        # get updated gradient (this will be normalized and passed in)
        grad = self.exp_avg / beta_correction1 / denom
        # renormalize and update warp
        # gradmax = self.eps + grad.reshape(grad.shape[0], -1).abs().max(1).values  # [B,]
        gradmax = self.eps + grad.norm(p=2, dim=-1, keepdim=True).flatten(1).max(1).values
        gradmax = gradmax.reshape(-1, *([1]) * (self.n_dims + 1))
        if not self.scaledown:  # if scaledown is "True", then we scale down even if the norm is below 1
            gradmax = torch.clamp(gradmax, min=1)
        # print(gradmax.abs().min(), gradmax.abs().max())
        grad = grad / gradmax * self.half_resolution  # norm is now 0.5r
        # multiply by learning rate
        grad.mul_(-self.lr)
        # print(grad.abs().max().item(), self.half_resolution, self.warp.shape)
        # compositional update
        w = grad + F.grid_sample(self.warp.data.permute(*self.permute_vtoimg), self.grid + grad, mode='bilinear',
                                 align_corners=True).permute(*self.permute_imgtov)
        # w = grad + self.warp.data
        # smooth result if asked for
        if self.smoothing_gaussians is not None:
            w = separable_filtering(w.permute(*self.permute_vtoimg), self.smoothing_gaussians).permute(
                *self.permute_imgtov)
        self.warp.data.copy_(w)
        # add to inverse if exists
        if self.optimize_inverse_warp and self.warpinv is not None:
            invwarp = compute_inverse_warp_displacement(self.warp.data, self.grid, self.warpinv.data, iters=5)
            warp_new = compute_inverse_warp_displacement(invwarp, self.grid, self.warp.data, iters=5)
            self.warp.data.copy_(warp_new)
            self.warpinv.data.copy_(invwarp)


class CompositiveWarp(nn.Module, AbstractDeformation):
    '''
    Class for compositive warp function (collects gradients of dL/dp)
    The image is computed as M \circ (\phi + u)
    '''

    def __init__(self, fixed_images: BatchedImages, moving_images: BatchedImages,
                 optimizer: str = 'Adam', optimizer_lr: float = 1e-2, optimizer_params: dict = {},
                 init_scale: int = 1,
                 smoothing_grad_sigma: float = 0.5, smoothing_warp_sigma: float = 0.5,
                 optimize_inverse_warp: bool = False,
                 ) -> None:
        super().__init__()
        self.num_images = num_images = fixed_images.size()
        spatial_dims = fixed_images.shape[2:]  # [H, W, [D]]
        self.n_dims = len(spatial_dims)  # number of spatial dimensions
        # permute indices
        self.permute_imgtov = (0, *range(2, self.n_dims + 2), 1)  # [N, HWD, dims] -> [N, HWD, dims] -> [N, dims, HWD]
        self.permute_vtoimg = (0, self.n_dims + 1, *range(1, self.n_dims + 1))  # [N, dims, HWD] -> [N, HWD, dims]
        self.device = fixed_images.device
        if optimizer_lr > 1:
            getLogger("CompositiveWarp").warning(
                f'optimizer_lr is {optimizer_lr}, which is very high. Unexpected registration may occur.')

        # define warp and register it as a parameter
        # define inverse warp and register it as a buffer
        self.optimize_inverse_warp = optimize_inverse_warp
        # set size
        if init_scale > 1:
            spatial_dims = [max(int(s / init_scale), MIN_IMG_SIZE) for s in spatial_dims]
        warp = torch.zeros([num_images, *spatial_dims, self.n_dims], dtype=torch.float32,
                           device=fixed_images.device)  # [N, HWD, dims]
        self.register_parameter('warp', nn.Parameter(warp))
        if self.optimize_inverse_warp:
            inv = torch.zeros([num_images, *spatial_dims, self.n_dims], dtype=torch.float32,
                              device=fixed_images.device)  # [N, HWD, dims]
        else:
            inv = torch.zeros([1], dtype=torch.float32, device=fixed_images.device)  # dummy
        self.register_buffer('inv', inv)

        # attach grad hook if smooothing of the gradient is required
        self.smoothing_grad_sigma = smoothing_grad_sigma
        if smoothing_grad_sigma > 0:
            self.smoothing_grad_gaussians = [gaussian_1d(s, truncated=2) for s in (
                        torch.zeros(self.n_dims, device=fixed_images.device) + smoothing_grad_sigma)]
        self.attach_grad_hook()

        # if the warp is also to be smoothed, add this constraint to the optimizer (in the optimizer_params dict)
        oparams = deepcopy(optimizer_params)
        self.smoothing_warp_sigma = smoothing_warp_sigma
        if self.smoothing_warp_sigma > 0:
            smoothing_warp_gaussians = [gaussian_1d(s, truncated=2) for s in
                                        (torch.zeros(self.n_dims, device=fixed_images.device) + smoothing_warp_sigma)]
            oparams['smoothing_gaussians'] = smoothing_warp_gaussians

        oparams['optimize_inverse_warp'] = optimize_inverse_warp
        if optimize_inverse_warp:
            oparams['warpinv'] = self.inv
        # add optimizer
        optimizer = optimizer.lower()
        if optimizer == 'sgd':
            raise Exception('lost by Eugenio when trimming down package')
            # self.optimizer = WarpSGD(self.warp, lr=optimizer_lr, **oparams)
        elif optimizer == 'adam':
            self.optimizer = WarpAdam(self.warp, lr=optimizer_lr, **oparams)
        else:
            raise NotImplementedError(f'Optimizer {optimizer} not implemented')

    def attach_grad_hook(self):
        ''' attack the grad hook to the velocity field if needed '''
        if self.smoothing_grad_sigma > 0:
            hook = partial(grad_smoothing_hook, gaussians=self.smoothing_grad_gaussians)
            self.warp.register_hook(hook)

    def initialize_grid(self):
        ''' initialize grid to a size
        Simply use the grid from the optimizer, which should be initialized to the correct size
        '''
        self.grid = self.optimizer.grid

    def set_zero_grad(self):
        ''' set the gradient to zero (or None) '''
        self.optimizer.zero_grad()

    def step(self):
        self.optimizer.step()

    def get_warp(self):
        ''' return warp function '''
        warp = self.warp
        return warp

    def get_inverse_warp(self, n_iters: int = 50, debug: bool = False, lr=0.1):
        ''' run an optimization procedure to get the inverse warp '''
        if self.optimize_inverse_warp:
            invwarp = self.inv
            invwarp = compute_inverse_warp_displacement(self.warp.data, self.grid, invwarp, iters=20)
        else:
            # no invwarp is defined, start from scratch
            invwarp = compute_inverse_warp_displacement(self.warp.data, self.grid, -self.warp.data, iters=200)
        return invwarp

    def set_size(self, size):
        # print(f"Setting size to {size}")
        ''' size: [H, W, D] or [H, W] '''
        mode = 'bilinear' if self.n_dims == 2 else 'trilinear'
        # get new displacement field
        warp = F.interpolate(self.warp.detach().permute(*self.permute_vtoimg), size=size, mode=mode, align_corners=True,
                             ).permute(*self.permute_imgtov)
        self.register_parameter('warp', nn.Parameter(warp))
        # set new inverse displacement field
        if len(self.inv.shape) > 1:
            self.inv = F.interpolate(self.inv.permute(*self.permute_vtoimg), size=size, mode=mode,
                                     align_corners=True).permute(*self.permute_imgtov)
        self.attach_grad_hook()
        self.optimizer.set_data_and_size(self.warp, size, warpinv=self.inv if self.optimize_inverse_warp else None)
        # interpolate inverse warp if it exists
        self.initialize_grid()


class GreedyRegistration(AbstractRegistration):
    '''
    This class implements greedy registration with a custom loss
    The moving image is interpolated to the fixed image grid, with an initial affine transform

    smooth_warp_sigma: how much to smooth the final warp field
    smooth_grad_sigma: how much to smooth the gradient of the final warp field  (this is similar to the Green's kernel)
    '''

    def __init__(self, scales: List[int], iterations: List[float],
                 fixed_images: BatchedImages, moving_images: BatchedImages,
                 loss_type: str = "cc",
                 deformation_type: str = 'geodesic',
                 optimizer: str = 'Adam', optimizer_params: dict = {},
                 optimizer_lr: float = 0.1,
                 integrator_n: Union[str, int] = 6,
                 mi_kernel_type: str = 'b-spline', cc_kernel_type: str = 'rectangular',
                 cc_kernel_size: int = 3,
                 smooth_warp_sigma: float = 0.5,
                 smooth_grad_sigma: float = 0.5,
                 loss_params: dict = {},
                 reduction: str = 'sum',
                 tolerance: float = 1e-6, max_tolerance_iters: int = 10,
                 init_affine: Optional[torch.Tensor] = None,
                 blur: bool = True,
                 custom_loss: nn.Module = None, **kwargs) -> None:
        # initialize abstract registration
        # nn.Module.__init__(self)
        super().__init__(scales=scales, iterations=iterations, fixed_images=fixed_images, moving_images=moving_images,
                         loss_type=loss_type, mi_kernel_type=mi_kernel_type, cc_kernel_type=cc_kernel_type,
                         custom_loss=custom_loss,
                         loss_params=loss_params,
                         cc_kernel_size=cc_kernel_size, reduction=reduction,
                         tolerance=tolerance, max_tolerance_iters=max_tolerance_iters, **kwargs)
        self.dims = fixed_images.dims
        self.blur = blur
        self.reduction = reduction
        # specify deformation type
        if deformation_type == 'geodesic':
            raise Exception('Eugenio lost this when trimming package')
            # warp = GeodesicShooting(fixed_images, moving_images, integrator_n=integrator_n, optimizer=optimizer,
            #                         optimizer_lr=optimizer_lr, optimizer_params=optimizer_params,
            #                         smoothing_grad_sigma=smooth_grad_sigma, init_scale=scales[0])
        elif deformation_type == 'compositive':
            warp = CompositiveWarp(fixed_images, moving_images, optimizer=optimizer, optimizer_lr=optimizer_lr,
                                   optimizer_params=optimizer_params, \
                                   smoothing_grad_sigma=smooth_grad_sigma, smoothing_warp_sigma=smooth_warp_sigma,
                                   init_scale=scales[0])
            smooth_warp_sigma = 0  # this work is delegated to compositive warp
        else:
            raise ValueError('Invalid deformation type: {}'.format(deformation_type))
        self.warp = warp
        self.smooth_warp_sigma = smooth_warp_sigma  # in voxels
        # initialize affine
        if init_affine is None:
            init_affine = torch.eye(self.dims + 1, device=fixed_images.device).unsqueeze(0).repeat(fixed_images.size(),
                                                                                                   1, 1)  # [N, D, D+1]
        self.affine = init_affine.detach()

    def get_warped_coordinates(self, fixed_images: BatchedImages, moving_images: BatchedImages, shape=None):
        ''' given fixed and moving images, get warp field (not displacement field) '''
        fixed_arrays = fixed_images()
        if shape is None:
            shape = fixed_images.shape
        else:
            shape = [fixed_arrays.shape[0], 1] + list(shape)

        fixed_t2p = fixed_images.get_torch2phy()
        moving_p2t = moving_images.get_phy2torch()
        # save initial affine transform to initialize grid
        affine_map_init = torch.matmul(moving_p2t, torch.matmul(self.affine, fixed_t2p))[:, :-1]
        # set affine coordinates
        fixed_image_affinecoords = F.affine_grid(affine_map_init, shape, align_corners=True)
        warp_field = self.warp.get_warp().clone()  # [N, HWD, 3]
        if tuple(warp_field.shape[1:-1]) != tuple(shape[2:]):
            # interpolate this
            warp_field = F.interpolate(warp_field.permute(*self.warp.permute_vtoimg), size=shape[2:], mode='trilinear',
                                       align_corners=True).permute(*self.warp.permute_imgtov)

        # smooth out the warp field if asked to
        if self.smooth_warp_sigma > 0:
            warp_gaussian = [gaussian_1d(s, truncated=2) for s in
                             (torch.zeros(self.dims, device=fixed_arrays.device) + self.smooth_warp_sigma)]
            warp_field = separable_filtering(warp_field.permute(*self.warp.permute_vtoimg), warp_gaussian).permute(
                *self.warp.permute_imgtov)
        # move these coordinates, and return them
        moved_coords = fixed_image_affinecoords + warp_field  # affine transform + warp field
        return moved_coords

    def evaluate(self, fixed_images: BatchedImages, moving_images: BatchedImages, shape=None):
        ''' given a new set of fixed and moving images, warp the fixed image '''
        moving_arrays = moving_images()
        moved_coords = self.get_warped_coordinates(fixed_images, moving_images, shape=shape)
        moved_image = F.grid_sample(moving_arrays, moved_coords, mode='bilinear',
                                    align_corners=True)  # [N, C, H, W, [D]]
        return moved_image

    def optimize(self, save_transformed=False):
        ''' optimize the warp field to match the two images based on loss function '''
        fixed_arrays = self.fixed_images()
        moving_arrays = self.moving_images()
        fixed_t2p = self.fixed_images.get_torch2phy()
        moving_p2t = self.moving_images.get_phy2torch()
        fixed_size = fixed_arrays.shape[2:]
        # save initial affine transform to initialize grid
        # init_grid = torch.eye(self.dims, self.dims+1).to(self.fixed_images.device).unsqueeze(0).repeat(self.fixed_images.size(), 1, 1)  # [N, dims, dims+1]
        affine_map_init = torch.matmul(moving_p2t, torch.matmul(self.affine, fixed_t2p))[:, :-1]

        # to save transformed images
        transformed_images = []
        # gaussian filter for smoothing the velocity field
        warp_gaussian = [gaussian_1d(s, truncated=2) for s in
                         (torch.zeros(self.dims, device=fixed_arrays.device) + self.smooth_warp_sigma)]
        # multi-scale optimization
        for scale, iters in zip(self.scales, self.iterations):
            self.convergence_monitor.reset()
            # resize images
            size_down = [max(int(s / scale), MIN_IMG_SIZE) for s in fixed_size]
            if self.blur and scale > 1:
                sigmas = 0.5 * torch.tensor([sz / szdown for sz, szdown in zip(fixed_size, size_down)],
                                            device=fixed_arrays.device)
                gaussians = [gaussian_1d(s, truncated=2) for s in sigmas]
                fixed_image_down = downsample(fixed_arrays, size=size_down, mode=self.fixed_images.interpolate_mode,
                                              gaussians=gaussians)
                moving_image_blur = separable_filtering(moving_arrays, gaussians)
            else:
                fixed_image_down = F.interpolate(fixed_arrays, size=size_down, mode=self.fixed_images.interpolate_mode,
                                                 align_corners=True)
                moving_image_blur = moving_arrays

            #### Set size for warp field
            self.warp.set_size(size_down)
            # Get coordinates to transform
            fixed_image_affinecoords = F.affine_grid(affine_map_init, fixed_image_down.shape, align_corners=True)
            pbar = tqdm(range(iters)) if self.progress_bar else range(iters)
            # reduce
            if self.reduction == 'mean':
                scale_factor = 1
            else:
                scale_factor = np.prod(fixed_image_down.shape)

            for i in pbar:
                self.warp.set_zero_grad()
                warp_field = self.warp.get_warp()  # [N, HWD, 3]
                # smooth out the warp field if asked to
                if self.smooth_warp_sigma > 0:
                    warp_field = separable_filtering(warp_field.permute(*self.warp.permute_vtoimg),
                                                     warp_gaussian).permute(*self.warp.permute_imgtov)
                moved_coords = fixed_image_affinecoords + warp_field  # affine transform + warp field
                # moved_coords.retain_grad()
                # move the image
                moved_image = F.grid_sample(moving_image_blur, moved_coords, mode='bilinear',
                                            align_corners=True)  # [N, C, H, W, [D]]
                loss = self.loss_fn(moved_image, fixed_image_down)
                loss.backward()
                if self.progress_bar:
                    pbar.set_description(
                        "scale: {}, iter: {}/{}, loss: {:4f}".format(scale, i, iters, loss.item() / scale_factor))
                # optimize the velocity field
                self.warp.step()
                # check for convergence
                if self.convergence_monitor.converged(loss.item()):
                    break

            # save transformed image
            if save_transformed:
                transformed_images.append(moved_image.detach().cpu())

        if save_transformed:
            return transformed_images


class HybridDiceLabelDiffloss(nn.Module):
    """
    Local squared zero-normalized cross-correlation.
    The loss is based on a moving kernel/window over the y_true/y_pred,
    within the window the square of zncc is calculated.
    The kernel can be a rectangular / triangular / gaussian window.
    The final loss is the averaged loss over all windows.
    Adapted from:
        https://github.com/voxelmorph/voxelmorph/blob/legacy/src/losses.py
        DeepReg (https://github.com/DeepRegNet/DeepReg)
    """

    def __init__(
            self,
            spatial_dims: int = 3,
            kernel_size: int = 3,
            kernel_type: str = "rectangular",
            reduction: str = "mean",
            smooth_nr: float = 1e-5,
            smooth_dr: float = 1e-5,
            unsigned: bool = True,
            checkpointing: bool = False,
            rel_weight_labeldiff: float = 1.0,
    ) -> None:
        """
        Args:
            spatial_dims: number of spatial dimensions, {``1``, ``2``, ``3``}. Defaults to 3.
            kernel_size: kernel spatial size, must be odd.
            kernel_type: {``"rectangular"``, ``"triangular"``, ``"gaussian"``}. Defaults to ``"rectangular"``.
            reduction: {``"none"``, ``"mean"``, ``"sum"``}
                Specifies the reduction to apply to the output. Defaults to ``"mean"``.
                - ``"none"``: no reduction will be applied.
                - ``"mean"``: the sum of the output will be divided by the number of elements in the output.
                - ``"sum"``: the output will be summed.
            smooth_nr: a small constant added to the numerator to avoid nan.
            smooth_dr: a small constant added to the denominator to avoid nan.
            rel_weight_labeldiff: relative weight of Label Difference loss
            split: do we want to split computation across 2 GPUs? (if pred and target are on different GPUs)
                default: False (assumes they are on same device and big enough to fit on one GPU)
        """
        super().__init__()
        self.ndim = spatial_dims
        if self.ndim != 3:
            raise ValueError("Unsupported ndim, only 3-d inputs are supported")
        if reduction != 'mean':
            raise ValueError("Unsupported reduction, only mean is supported")
        self.reduction = reduction
        self.unsigned = unsigned

        self.kernel_size = kernel_size
        if self.kernel_size % 2 == 0:
            raise ValueError(f"kernel_size must be odd, got {self.kernel_size}")

        # _kernel = look_up_option(kernel_type, kernel_dict)
        _kernel = kernel_dict[kernel_type]
        self.kernel = _kernel(self.kernel_size)
        self.kernel.requires_grad = False
        self.kernel_nd, self.kernel_vol = self.get_kernel_vol()  # get nD kernel and its volume
        self.smooth_nr = float(smooth_nr)
        self.smooth_dr = float(smooth_dr)
        self.rel_weight_labeldiff = rel_weight_labeldiff
        self.checkpointing = checkpointing

    def get_kernel_vol(self):
        vol = self.kernel
        for _ in range(self.ndim - 1):
            vol = torch.matmul(vol.unsqueeze(-1), self.kernel.unsqueeze(0))
        return vol, torch.sum(vol)

    def forward(self, pred: torch.Tensor, target: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Args:
            pred: the shape should be BNH[WD].
            target: the shape should be BNH[WD].
        Raises:
            ValueError: When ``self.reduction`` is not "mean"".
        """
        if pred.ndim - 2 != 3:
            raise ValueError(f"expecting pred with 3 spatial dimensions, got pred of shape {pred.shape}")
        if target.shape != pred.shape:
            raise ValueError(f"ground truth has differing shape ({target.shape}) from pred ({pred.shape})")
        if mask is not None:
            raise ValueError('Mask should always be None')

        # sum over kernel
        def cc_checkpoint_fn(target, pred, kernel, kernel_vol):
            '''
            This function is used to compute the intermediate results of the loss.
            '''
            t2, p2, tp = target * target, pred * pred, target * pred
            kernel, kernel_vol = kernel.to(pred), kernel_vol.to(pred)
            # kernel_nd = self.kernel_nd.to(pred)
            kernels = [kernel] * self.ndim
            kernels_t = kernels_p = kernels
            kernel_vol_t = kernel_vol_p = kernel_vol
            # compute intermediates
            t_sum = separable_filtering(target, kernels=kernels_t)
            p_sum = separable_filtering(pred, kernels=kernels_p)
            t2_sum = separable_filtering(t2, kernels=kernels_t)
            p2_sum = separable_filtering(p2, kernels=kernels_p)
            tp_sum = separable_filtering(tp, kernels=kernels_t)  # use target device's output
            # average over kernel
            t_avg = t_sum / kernel_vol_t
            p_avg = p_sum / kernel_vol_p
            # normalized cross correlation between t and p
            # sum[(t - mean[t]) * (p - mean[p])] / std[t] / std[p]
            # denoted by num / denom
            # assume we sum over N values
            # num = sum[t * p - mean[t] * p - t * mean[p] + mean[t] * mean[p]]
            #     = sum[t*p] - sum[t] * sum[p] / N * 2 + sum[t] * sum[p] / N
            #     = sum[t*p] - sum[t] * sum[p] / N
            #     = sum[t*p] - sum[t] * mean[p] = cross
            # the following is actually squared ncc
            cross = (tp_sum.to(pred) - p_avg * t_sum.to(pred))  # on pred device
            t_var = torch.max(
                t2_sum - t_avg * t_sum, torch.as_tensor(self.smooth_dr, dtype=t2_sum.dtype, device=t2_sum.device)
            ).to(pred)
            p_var = torch.max(
                p2_sum - p_avg * p_sum, torch.as_tensor(self.smooth_dr, dtype=p2_sum.dtype, device=p2_sum.device)
            )
            if self.unsigned:
                ncc: torch.Tensor = (cross * cross + self.smooth_nr) / ((t_var * p_var) + self.smooth_dr)
            else:
                ncc: torch.Tensor = (cross + self.smooth_nr) / (
                            (torch.sqrt(t_var) * torch.sqrt(p_var)) + self.smooth_dr)
            return ncc

        if self.checkpointing:
            ncc = checkpoint(cc_checkpoint_fn, target[:, 0:1, :, :, :], pred[:, 0:1, :, :, :], self.kernel,
                             self.kernel_vol)
        else:
            ncc = cc_checkpoint_fn(target[:, 0:1, :, :, :], pred[:, 0:1, :, :, :], self.kernel, self.kernel_vol)

        lnccloss = (1.0 - torch.mean(ncc))

        diff = torch.abs(target[:, 1:target.shape[1], :, :, :] - pred[:, 1:target.shape[1], :, :, :])
        labeldiffloss = 0.5 * diff.sum(dim=1).mean([1, 2, 3])

        loss = lnccloss + labeldiffloss * self.rel_weight_labeldiff

        return (loss / (1.0 + self.rel_weight_labeldiff))
