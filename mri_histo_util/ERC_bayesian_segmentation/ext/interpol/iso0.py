"""Isotropic 0-th order splines ("nearest neighbor")"""
import torch
from .bounds import Bound
from .jit_utils import (sub2ind_list, make_sign,
                        inbounds_mask_3d, inbounds_mask_2d, inbounds_mask_1d)
from typing import List, Optional
Tensor = torch.Tensor


@torch.jit.script
def get_indices(g, n: int, bound: Bound):
    g0 = g.round().long()
    sign0 = bound.transform(g0, n)
    g0 = bound.index(g0, n)
    return g0, sign0


# ======================================================================
#                                 3D
# ======================================================================


@torch.jit.script
def pull3d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, iZ) tensor
    g: (B, oX, oY, oZ, 3) tensor
    bound: List{3}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY, oZ) tensor
    """
    dim = 3
    boundx, boundy, boundz = bound
    oshape = g.shape[-dim-1:-1]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = g.unbind(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = inp.shape[-dim:]
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)
    gy, signy = get_indices(gy, ny, boundy)
    gz, signz = get_indices(gz, nz, boundz)

    # gather
    inp = inp.reshape(inp.shape[:2] + [-1])
    idx = sub2ind_list([gx, gy, gz], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = make_sign([signx, signy, signz])
    if sign is not None:
        out *= sign
    if mask is not None:
        out *= mask
    out = out.reshape(out.shape[:2] + oshape)
    return out


@torch.jit.script
def push3d(inp, g, shape: Optional[List[int]], bound: List[Bound],
           extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, iZ) tensor
    g: (B, iX, iY, iZ, 3) tensor
    shape: List{3}[int], optional
    bound: List{3}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 3
    boundx, boundy, boundz = bound
    if inp.shape[-dim:] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = inp.shape[-dim:]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = torch.unbind(g, -1)
    inp = inp.reshape(inp.shape[:2] + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)
    gy, signy = get_indices(gy, ny, boundy)
    gz, signz = get_indices(gz, nz, boundz)

    # scatter
    out = torch.zeros([batch, channel, nx*ny*nz], dtype=inp.dtype, device=inp.device)
    idx = sub2ind_list([gx, gy, gz], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    sign = make_sign([signx, signy, signz])
    if sign is not None or mask is not None:
        inp = inp.clone()
    if sign is not None:
        inp *= sign
    if mask is not None:
        inp *= mask
    out.scatter_add_(-1, idx, inp)

    out = out.reshape(out.shape[:2] + shape)
    return out


# ======================================================================
#                                 2D
# ======================================================================


@torch.jit.script
def pull2d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY) tensor
    g: (B, oX, oY, 2) tensor
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY) tensor
    """
    dim = 2
    boundx, boundy = bound
    oshape = g.shape[-dim-1:-1]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = g.unbind(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = inp.shape[-dim:]
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)
    gy, signy = get_indices(gy, ny, boundy)

    # gather
    inp = inp.reshape(inp.shape[:2] + [-1])
    idx = sub2ind_list([gx, gy], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = make_sign([signx, signy])
    if sign is not None:
        out = out * sign
    if mask is not None:
        out = mask * mask
    out = out.reshape(out.shape[:2] + oshape)
    return out


@torch.jit.script
def push2d(inp, g, shape: Optional[List[int]], bound: List[Bound],
           extrapolate: int = 1):
    """
    inp: (B, C, iX, iY) tensor
    g: (B, iX, iY, 2) tensor
    shape: List{2}[int], optional
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 2
    boundx, boundy = bound
    if inp.shape[-dim:] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = inp.shape[-dim:]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = torch.unbind(g, -1)
    inp = inp.reshape(inp.shape[:2] + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)
    gy, signy = get_indices(gy, ny, boundy)

    # scatter
    out = torch.zeros([batch, channel, nx*ny], dtype=inp.dtype, device=inp.device)
    idx = sub2ind_list([gx, gy], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    sign = make_sign([signx, signy])
    if sign is not None or mask is not None:
        inp = inp.clone()
    if sign is not None:
        inp = inp * sign
    if mask is not None:
        inp = inp * mask
    out.scatter_add_(-1, idx, inp)

    out = out.reshape(out.shape[:2] + shape)
    return out


# ======================================================================
#                                 1D
# ======================================================================


@torch.jit.script
def pull1d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX) tensor
    g: (B, oX, 1) tensor
    bound: List{1}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX) tensor
    """
    dim = 1
    boundx = bound[0]
    oshape = g.shape[-dim-1:-1]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = inp.shape[-dim:]
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)

    # gather
    inp = inp.reshape(inp.shape[:2] + [-1])
    idx = gx
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = signx
    if sign is not None:
        out = out * sign
    if mask is not None:
        out = out * mask
    out = out.reshape(out.shape[:2] + oshape)
    return out


@torch.jit.script
def push1d(inp, g, shape: Optional[List[int]], bound: List[Bound],
           extrapolate: int = 1):
    """
    inp: (B, C, iX) tensor
    g: (B, iX, 1) tensor
    shape: List{1}[int], optional
    bound: List{1}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 1
    boundx = bound[0]
    if inp.shape[-dim:] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = inp.shape[-dim:]
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    inp = inp.reshape(inp.shape[:2] + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # nearest integer coordinates
    gx, signx = get_indices(gx, nx, boundx)

    # scatter
    out = torch.zeros([batch, channel, nx], dtype=inp.dtype, device=inp.device)
    idx = gx
    idx = idx.expand([batch, channel, idx.shape[-1]])
    sign = signx
    if sign is not None or mask is not None:
        inp = inp.clone()
    if sign is not None:
        inp = inp * sign
    if mask is not None:
        inp = inp * mask
    out.scatter_add_(-1, idx, inp)

    out = out.reshape(out.shape[:2] + shape)
    return out


# ======================================================================
#                                 ND
# ======================================================================


@torch.jit.script
def grad(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, *ishape) tensor
    g: (B, *oshape, D) tensor
    bound: List{D}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *oshape, D) tensor
    """
    dim = g.shape[-1]
    oshape = list(g.shape[-dim-1:-1])
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    return torch.zeros([batch, channel] + oshape + [dim],
                       dtype=inp.dtype, device=inp.device)


@torch.jit.script
def pushgrad(inp, g, shape: Optional[List[int]], bound: List[Bound],
             extrapolate: int = 1):
    """
    inp: (B, C, *ishape, D) tensor
    g: (B, *ishape, D) tensor
    shape: List{D}[int], optional, optional
    bound: List{D}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = g.shape[-1]
    if inp.shape[-dim-1:-1] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = inp.shape[-dim-1:-1]
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)

    return torch.zeros([batch, channel] + shape,
                       dtype=inp.dtype, device=inp.device)


@torch.jit.script
def hess(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, *ishape) tensor
    g: (B, *oshape, D) tensor
    bound: List{D}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *oshape, D, D) tensor
    """
    dim = g.shape[-1]
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    return torch.zeros([batch, channel] + oshape + [dim, dim],
                       dtype=inp.dtype, device=inp.device)
