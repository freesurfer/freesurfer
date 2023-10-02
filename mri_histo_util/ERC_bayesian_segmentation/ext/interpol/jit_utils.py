"""A lot of utility functions for TorchScript"""
import torch
import os
from typing import List, Tuple, Optional
from .utils import torch_version
from torch import Tensor


@torch.jit.script
def pad_list_int(x: List[int], dim: int) -> List[int]:
    if len(x) < dim:
        x = x + x[-1:] * (dim - len(x))
    if len(x) > dim:
        x = x[:dim]
    return x


@torch.jit.script
def pad_list_float(x: List[float], dim: int) -> List[float]:
    if len(x) < dim:
        x = x + x[-1:] * (dim - len(x))
    if len(x) > dim:
        x = x[:dim]
    return x


@torch.jit.script
def pad_list_str(x: List[str], dim: int) -> List[str]:
    if len(x) < dim:
        x = x + x[-1:] * (dim - len(x))
    if len(x) > dim:
        x = x[:dim]
    return x


@torch.jit.script
def list_any(x: List[bool]) -> bool:
    for elem in x:
        if elem:
            return True
    return False


@torch.jit.script
def list_all(x: List[bool]) -> bool:
    for elem in x:
        if not elem:
            return False
    return True


@torch.jit.script
def list_prod_int(x: List[int]) -> int:
    if len(x) == 0:
        return 1
    x0 = x[0]
    for x1 in x[1:]:
        x0 = x0 * x1
    return x0


@torch.jit.script
def list_sum_int(x: List[int]) -> int:
    if len(x) == 0:
        return 1
    x0 = x[0]
    for x1 in x[1:]:
        x0 = x0 + x1
    return x0


@torch.jit.script
def list_prod_tensor(x: List[Tensor]) -> Tensor:
    if len(x) == 0:
        empty: List[int] = []
        return torch.ones(empty)
    x0 = x[0]
    for x1 in x[1:]:
        x0 = x0 * x1
    return x0


@torch.jit.script
def list_sum_tensor(x: List[Tensor]) -> Tensor:
    if len(x) == 0:
        empty: List[int] = []
        return torch.ones(empty)
    x0 = x[0]
    for x1 in x[1:]:
        x0 = x0 + x1
    return x0


@torch.jit.script
def list_reverse_int(x: List[int]) -> List[int]:
    if len(x) == 0:
        return x
    return [x[i] for i in range(-1, -len(x)-1, -1)]


@torch.jit.script
def list_cumprod_int(x: List[int], reverse: bool = False,
                     exclusive: bool = False) -> List[int]:
    if len(x) == 0:
        lx: List[int] = []
        return lx
    if reverse:
        x = list_reverse_int(x)

    x0 = 1 if exclusive else x[0]
    lx = [x0]
    all_x = x[:-1] if exclusive else x[1:]
    for x1 in all_x:
        x0 = x0 * x1
        lx.append(x0)
    if reverse:
        lx = list_reverse_int(lx)
    return lx


@torch.jit.script
def movedim1(x, source: int, destination: int):
    dim = x.dim()
    source = dim + source if source < 0 else source
    destination = dim + destination if destination < 0 else destination
    permutation = [d for d in range(dim)]
    permutation = permutation[:source] + permutation[source+1:]
    permutation = permutation[:destination] + [source] + permutation[destination:]
    return x.permute(permutation)


@torch.jit.script
def sub2ind(subs, shape: List[int]):
    """Convert sub indices (i, j, k) into linear indices.

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    subs : (D, ...) tensor
        List of sub-indices. The first dimension is the number of dimension.
        Each element should have the same number of elements and shape.
    shape : (D,) list[int]
        Size of each dimension. Its length should be the same as the
        first dimension of ``subs``.

    Returns
    -------
    ind : (...) tensor
        Linear indices
    """
    subs = subs.unbind(0)
    ind = subs[-1]
    subs = subs[:-1]
    ind = ind.clone()
    stride = list_cumprod_int(shape[1:], reverse=True, exclusive=False)
    for i, s in zip(subs, stride):
        ind += i * s
    return ind


@torch.jit.script
def sub2ind_list(subs: List[Tensor], shape: List[int]):
    """Convert sub indices (i, j, k) into linear indices.

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    subs : (D,) list[tensor]
        List of sub-indices. The first dimension is the number of dimension.
        Each element should have the same number of elements and shape.
    shape : (D,) list[int]
        Size of each dimension. Its length should be the same as the
        first dimension of ``subs``.

    Returns
    -------
    ind : (...) tensor
        Linear indices
    """
    ind = subs[-1]
    subs = subs[:-1]
    ind = ind.clone()
    stride = list_cumprod_int(shape[1:], reverse=True, exclusive=False)
    for i, s in zip(subs, stride):
        ind += i * s
    return ind

# floor_divide returns wrong results for negative values, because it truncates
# instead of performing a proper floor. In recent version of pytorch, it is
# advised to use div(..., rounding_mode='trunc'|'floor') instead.
# Here, we only use floor_divide on positive values so we do not care.
if torch_version('>=', [1, 8]):
    @torch.jit.script
    def floor_div(x, y) -> torch.Tensor:
        return torch.div(x, y, rounding_mode='floor')
    @torch.jit.script
    def floor_div_int(x, y: int) -> torch.Tensor:
        return torch.div(x, y, rounding_mode='floor')
else:
    @torch.jit.script
    def floor_div(x, y) -> torch.Tensor:
        return (x / y).floor_()
    @torch.jit.script
    def floor_div_int(x, y: int) -> torch.Tensor:
        return (x / y).floor_()


@torch.jit.script
def ind2sub(ind, shape: List[int]):
    """Convert linear indices into sub indices (i, j, k).

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    ind : tensor_like
        Linear indices
    shape : (D,) vector_like
        Size of each dimension.

    Returns
    -------
    subs : (D, ...) tensor
        Sub-indices.
    """
    stride = list_cumprod_int(shape, reverse=True, exclusive=True)
    sub = ind.new_empty([len(shape)] + ind.shape)
    sub.copy_(ind)
    for d in range(len(shape)):
        if d > 0:
            sub[d] = torch.remainder(sub[d], stride[d-1])
        sub[d] = floor_div_int(sub[d], stride[d])
    return sub


@torch.jit.script
def inbounds_mask_3d(extrapolate: int, gx, gy, gz, nx: int, ny: int, nz: int) \
        -> Optional[Tensor]:
    # mask of inbounds voxels
    mask: Optional[Tensor] = None
    if extrapolate in (0, 2):  # no / hist
        tiny = 5e-2
        threshold = tiny
        if extrapolate == 2:
            threshold = 0.5 + tiny
        mask = ((gx > -threshold) & (gx < nx - 1 + threshold) &
                (gy > -threshold) & (gy < ny - 1 + threshold) &
                (gz > -threshold) & (gz < nz - 1 + threshold))
        return mask
    return mask


@torch.jit.script
def inbounds_mask_2d(extrapolate: int, gx, gy, nx: int, ny: int) \
        -> Optional[Tensor]:
    # mask of inbounds voxels
    mask: Optional[Tensor] = None
    if extrapolate in (0, 2):  # no / hist
        tiny = 5e-2
        threshold = tiny
        if extrapolate == 2:
            threshold = 0.5 + tiny
        mask = ((gx > -threshold) & (gx < nx - 1 + threshold) &
                (gy > -threshold) & (gy < ny - 1 + threshold))
        return mask
    return mask


@torch.jit.script
def inbounds_mask_1d(extrapolate: int, gx, nx: int) -> Optional[Tensor]:
    # mask of inbounds voxels
    mask: Optional[Tensor] = None
    if extrapolate in (0, 2):  # no / hist
        tiny = 5e-2
        threshold = tiny
        if extrapolate == 2:
            threshold = 0.5 + tiny
        mask = (gx > -threshold) & (gx < nx - 1 + threshold)
        return mask
    return mask


@torch.jit.script
def make_sign(sign: List[Optional[Tensor]]) -> Optional[Tensor]:
    is_none : List[bool] = [s is None for s in sign]
    if list_all(is_none):
        return None
    filt_sign: List[Tensor] = []
    for s in sign:
        if s is not None:
            filt_sign.append(s)
    return list_prod_tensor(filt_sign)


@torch.jit.script
def square(x):
    return x * x


@torch.jit.script
def square_(x):
    return x.mul_(x)


@torch.jit.script
def cube(x):
    return x * x * x


@torch.jit.script
def cube_(x):
    return square_(x).mul_(x)


@torch.jit.script
def pow4(x):
    return square(square(x))


@torch.jit.script
def pow4_(x):
    return square_(square_(x))


@torch.jit.script
def pow5(x):
    return x * pow4(x)


@torch.jit.script
def pow5_(x):
    return pow4_(x).mul_(x)


@torch.jit.script
def pow6(x):
    return square(cube(x))


@torch.jit.script
def pow6_(x):
    return square_(cube_(x))


@torch.jit.script
def pow7(x):
    return pow6(x) * x


@torch.jit.script
def pow7_(x):
    return pow6_(x).mul_(x)


@torch.jit.script
def dot(x, y, dim: int = -1, keepdim: bool = False):
    """(Batched) dot product along a dimension"""
    x = movedim1(x, dim, -1).unsqueeze(-2)
    y = movedim1(y, dim, -1).unsqueeze(-1)
    d = torch.matmul(x, y).squeeze(-1).squeeze(-1)
    if keepdim:
        d.unsqueeze(dim)
    return d


@torch.jit.script
def dot_multi(x, y, dim: List[int], keepdim: bool = False):
    """(Batched) dot product along a dimension"""
    for d in dim:
        x = movedim1(x, d, -1)
        y = movedim1(y, d, -1)
    x = x.reshape(x.shape[:-len(dim)] + [1, -1])
    y = y.reshape(x.shape[:-len(dim)] + [-1, 1])
    dt = torch.matmul(x, y).squeeze(-1).squeeze(-1)
    if keepdim:
        for d in dim:
            dt.unsqueeze(d)
    return dt



# cartesian_prod takes multiple inout tensors as input in eager mode
# but takes a list of tensor in jit mode. This is a helper that works
# in both cases.
if not int(os.environ.get('PYTORCH_JIT', '1')):
    cartesian_prod = lambda x: torch.cartesian_prod(*x)
    if torch_version('>=', (1, 10)):
        def meshgrid_ij(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(*x, indexing='ij')
        def meshgrid_xy(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(*x, indexing='xy')
    else:
        def meshgrid_ij(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(*x)
        def meshgrid_xy(x: List[torch.Tensor]) -> List[torch.Tensor]:
            grid = torch.meshgrid(*x)
            if len(grid) > 1:
                grid[0] = grid[0].transpose(0, 1)
                grid[1] = grid[1].transpose(0, 1)
            return grid

else:
    cartesian_prod = torch.cartesian_prod
    if torch_version('>=', (1, 10)):
        @torch.jit.script
        def meshgrid_ij(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(x, indexing='ij')
        @torch.jit.script
        def meshgrid_xy(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(x, indexing='xy')
    else:
        @torch.jit.script
        def meshgrid_ij(x: List[torch.Tensor]) -> List[torch.Tensor]:
            return torch.meshgrid(x)
        @torch.jit.script
        def meshgrid_xyt(x: List[torch.Tensor]) -> List[torch.Tensor]:
            grid = torch.meshgrid(x)
            if len(grid) > 1:
                grid[0] = grid[0].transpose(0, 1)
                grid[1] = grid[1].transpose(0, 1)
            return grid


meshgrid = meshgrid_ij


# In torch < 1.6, div applied to integer tensor performed a floor_divide
# In torch > 1.6, it performs a true divide.
# Floor division must be done using `floor_divide`, but it was buggy
# until torch 1.13 (it was doing a trunc divide instead of a floor divide).
# There was at some point a deprecation warning for floor_divide, but it
# seems to have been lifted afterwards. In torch >= 1.13, floor_divide
# performs a correct floor division.
# Since we only apply floor_divide ot positive values, we are fine.
if torch_version('<', (1, 6)):
    floor_div = torch.div
else:
    floor_div = torch.floor_divide