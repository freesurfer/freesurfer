"""
Non-differentiable forward/backward components.
These components are put together in `interpol.autograd` to generate
differentiable functions.

Note
----
.. I removed @torch.jit.script from these entry-points because compiling
   all possible combinations of bound+interpolation made the first call
   extremely slow.
.. I am not using the dot/multi_dot helpers even though they should be
   more efficient that "multiply and sum" because I haven't had the time
   to test them. It would be worth doing it.
"""
import torch
from typing import List, Optional, Tuple
from .jit_utils import list_all, dot, dot_multi, pad_list_int
from .bounds import Bound
from .splines import Spline
from . import iso0, iso1, nd
Tensor = torch.Tensor


@torch.jit.script
def make_bound(bound: List[int]) -> List[Bound]:
    return [Bound(b) for b in bound]


@torch.jit.script
def make_spline(spline: List[int]) -> List[Spline]:
    return [Spline(s) for s in spline]


# @torch.jit.script
def grid_pull(inp, grid, bound: List[int], interpolation: List[int],
              extrapolate: int):
    """
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_out, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_out) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.pull3d(inp, grid, bound_fn, extrapolate)
        elif dim == 2:
            return iso1.pull2d(inp, grid, bound_fn, extrapolate)
        elif dim == 1:
            return iso1.pull1d(inp, grid, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        if dim == 3:
            return iso0.pull3d(inp, grid, bound_fn, extrapolate)
        elif dim == 2:
            return iso0.pull2d(inp, grid, bound_fn, extrapolate)
        elif dim == 1:
            return iso0.pull1d(inp, grid, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.pull(inp, grid, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_push(inp, grid, shape: Optional[List[int]], bound: List[int],
              interpolation: List[int], extrapolate: int):
    """
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_in, D) tensor
    shape: List{D}[int] tensor, optional, default=spatial_in
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *shape) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.push3d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 2:
            return iso1.push2d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 1:
            return iso1.push1d(inp, grid, shape, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        if dim == 3:
            return iso0.push3d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 2:
            return iso0.push2d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 1:
            return iso0.push1d(inp, grid, shape, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.push(inp, grid, shape, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_count(grid, shape: Optional[List[int]], bound: List[int],
              interpolation: List[int], extrapolate: int):
    """
    grid: (B, *spatial_in, D) tensor
    shape: List{D}[int] tensor, optional, default=spatial_in
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, 1, *shape) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    gshape = list(grid.shape[-dim-1:-1])
    if shape is None:
        shape = gshape
    inp = torch.ones([], dtype=grid.dtype, device=grid.device)
    inp = inp.expand([len(grid), 1] + gshape)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.push3d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 2:
            return iso1.push2d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 1:
            return iso1.push1d(inp, grid, shape, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        if dim == 3:
            return iso0.push3d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 2:
            return iso0.push2d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 1:
            return iso0.push1d(inp, grid, shape, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.push(inp, grid, shape, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_grad(inp, grid, bound: List[int], interpolation: List[int],
              extrapolate: int):
    """
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_out, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_out, D) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.grad3d(inp, grid, bound_fn, extrapolate)
        elif dim == 2:
            return iso1.grad2d(inp, grid, bound_fn, extrapolate)
        elif dim == 1:
            return iso1.grad1d(inp, grid, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        return iso0.grad(inp, grid, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.grad(inp, grid, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_pushgrad(inp, grid, shape: List[int], bound: List[int],
                  interpolation: List[int], extrapolate: int):
    """ /!\ Used only in backward pass of grid_grad
    inp: (B, C, *spatial_in, D) tensor
    grid: (B, *spatial_in, D) tensor
    shape: List{D}[int], optional
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *shape) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.pushgrad3d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 2:
            return iso1.pushgrad2d(inp, grid, shape, bound_fn, extrapolate)
        elif dim == 1:
            return iso1.pushgrad1d(inp, grid, shape, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        return iso0.pushgrad(inp, grid, shape, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.pushgrad(inp, grid, shape, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_hess(inp, grid, bound: List[int], interpolation: List[int],
              extrapolate: int):
    """ /!\ Used only in backward pass of grid_grad
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_out, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_out, D, D) tensor
    """
    dim = grid.shape[-1]
    bound = pad_list_int(bound, dim)
    interpolation = pad_list_int(interpolation, dim)
    bound_fn = make_bound(bound)
    is_iso1 = list_all([order == 1 for order in interpolation])
    if is_iso1:
        if dim == 3:
            return iso1.hess3d(inp, grid, bound_fn, extrapolate)
        if dim == 2:
            return iso1.hess2d(inp, grid, bound_fn, extrapolate)
        if dim == 1:
            return iso1.hess1d(inp, grid, bound_fn, extrapolate)
    is_iso0 = list_all([order == 0 for order in interpolation])
    if is_iso0:
        return iso0.hess(inp, grid, bound_fn, extrapolate)
    spline_fn = make_spline(interpolation)
    return nd.hess(inp, grid, bound_fn, spline_fn, extrapolate)


# @torch.jit.script
def grid_pull_backward(grad, inp, grid, bound: List[int],
                       interpolation: List[int], extrapolate: int) \
        -> Tuple[Optional[Tensor], Optional[Tensor], ]:
    """
    grad: (B, C, *spatial_out) tensor
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_out, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_in) tensor, (B, *spatial_out, D)
    """
    dim = grid.shape[-1]
    grad_inp: Optional[Tensor] = None
    grad_grid: Optional[Tensor] = None
    if inp.requires_grad:
        grad_inp = grid_push(grad, grid, inp.shape[-dim:], bound, interpolation, extrapolate)
    if grid.requires_grad:
        grad_grid = grid_grad(inp, grid, bound, interpolation, extrapolate)
        # grad_grid = dot(grad_grid, grad.unsqueeze(-1), dim=1)
        grad_grid = (grad_grid * grad.unsqueeze(-1)).sum(dim=1)
    return grad_inp, grad_grid


# @torch.jit.script
def grid_push_backward(grad, inp, grid, bound: List[int],
                       interpolation: List[int], extrapolate: int) \
        -> Tuple[Optional[Tensor], Optional[Tensor], ]:
    """
    grad: (B, C, *spatial_out) tensor
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_in, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_in) tensor, (B, *spatial_in, D)
    """
    grad_inp: Optional[Tensor] = None
    grad_grid: Optional[Tensor] = None
    if inp.requires_grad:
        grad_inp = grid_pull(grad, grid, bound, interpolation, extrapolate)
    if grid.requires_grad:
        grad_grid = grid_grad(grad, grid, bound, interpolation, extrapolate)
        # grad_grid = dot(grad_grid, inp.unsqueeze(-1), dim=1)
        grad_grid = (grad_grid * inp.unsqueeze(-1)).sum(dim=1)
    return grad_inp, grad_grid


# @torch.jit.script
def grid_count_backward(grad, grid, bound: List[int],
                       interpolation: List[int], extrapolate: int) \
        -> Optional[Tensor]:
    """
    grad: (B, C, *spatial_out) tensor
    grid: (B, *spatial_in, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_in) tensor, (B, *spatial_in, D)
    """
    if grid.requires_grad:
        return grid_grad(grad, grid, bound, interpolation, extrapolate).sum(1)
    return None


# @torch.jit.script
def grid_grad_backward(grad, inp, grid, bound: List[int],
                       interpolation: List[int], extrapolate: int) \
        -> Tuple[Optional[Tensor], Optional[Tensor]]:
    """
    grad: (B, C, *spatial_out, D) tensor
    inp: (B, C, *spatial_in) tensor
    grid: (B, *spatial_out, D) tensor
    bound: List{D}[int] tensor
    interpolation: List{D}[int]
    extrapolate: int
    returns: (B, C, *spatial_in, D) tensor, (B, *spatial_out, D)
    """
    dim = grid.shape[-1]
    shape = inp.shape[-dim:]
    grad_inp: Optional[Tensor] = None
    grad_grid: Optional[Tensor] = None
    if inp.requires_grad:
        grad_inp = grid_pushgrad(grad, grid, shape, bound, interpolation, extrapolate)
    if grid.requires_grad:
        grad_grid = grid_hess(inp, grid, bound, interpolation, extrapolate)
        # grad_grid = dot_multi(grad_grid, grad.unsqueeze(-1), dim=[1, -2])
        grad_grid = (grad_grid * grad.unsqueeze(-1)).sum(dim=[1, -2])
    return grad_inp, grad_grid
