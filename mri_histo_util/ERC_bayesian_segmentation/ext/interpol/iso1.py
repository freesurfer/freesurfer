"""Isotropic 1-st order splines ("linear/bilinear/trilinear")"""
import torch
from .bounds import Bound
from .jit_utils import (sub2ind_list, make_sign,
                        inbounds_mask_3d, inbounds_mask_2d, inbounds_mask_1d)
from typing import List, Tuple, Optional
Tensor = torch.Tensor


@torch.jit.script
def get_weights_and_indices(g, n: int, bound: Bound) \
        -> Tuple[Tensor, Tensor, Tensor, Optional[Tensor], Optional[Tensor]]:
    g0 = g.floor().long()
    g1 = g0 + 1
    sign1 = bound.transform(g1, n)
    sign0 = bound.transform(g0, n)
    g1 = bound.index(g1, n)
    g0 = bound.index(g0, n)
    g = g - g.floor()
    return g, g0, g1, sign0, sign1


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
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = g.unbind(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)
    gz, gz0, gz1, signz0, signz1 = get_weights_and_indices(gz, nz, boundz)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    # - corner 000
    idx = sub2ind_list([gx0, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = make_sign([signx0, signy0, signz0])
    if sign is not None:
        out = out * sign
    out = out * ((1 - gx) * (1 - gy) * (1 - gz))
    # - corner 001
    idx = sub2ind_list([gx0, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy0, signz1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * ((1 - gx) * (1 - gy) * gz)
    out = out + out1
    # - corner 010
    idx = sub2ind_list([gx0, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz0])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * ((1 - gx) * gy * (1 - gz))
    out = out + out1
    # - corner 011
    idx = sub2ind_list([gx0, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * ((1 - gx) * gy * gz)
    out = out + out1
    # - corner 100
    idx = sub2ind_list([gx1, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz0])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * (1 - gy) * (1 - gz))
    out = out + out1
    # - corner 101
    idx = sub2ind_list([gx1, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * (1 - gy) * gz)
    out = out + out1
    # - corner 110
    idx = sub2ind_list([gx1, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz0])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * gy * (1 - gz))
    out = out + out1
    # - corner 111
    idx = sub2ind_list([gx1, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * gy * gz)
    out = out + out1

    if mask is not None:
        out *= mask
    out = out.reshape(list(out.shape[:2]) + oshape)
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
    ishape = list(inp.shape[-dim:])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = torch.unbind(g, -1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)
    gz, gz0, gz1, signz0, signz1 = get_weights_and_indices(gz, nz, boundz)

    # scatter
    out = torch.zeros([batch, channel, nx*ny*nz],
                      dtype=inp.dtype, device=inp.device)
    # - corner 000
    idx = sub2ind_list([gx0, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0, signz0])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * ((1 - gx) * (1 - gy) * (1 - gz))
    out.scatter_add_(-1, idx, out1)
    # - corner 001
    idx = sub2ind_list([gx0, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0, signz1])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * ((1 - gx) * (1 - gy) * gz)
    out.scatter_add_(-1, idx, out1)
    # - corner 010
    idx = sub2ind_list([gx0, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1, signz0])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * ((1 - gx) * gy * (1 - gz))
    out.scatter_add_(-1, idx, out1)
    # - corner 011
    idx = sub2ind_list([gx0, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1, signz1])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * ((1 - gx) * gy * gz)
    out.scatter_add_(-1, idx, out1)
    # - corner 100
    idx = sub2ind_list([gx1, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0, signz0])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * (gx * (1 - gy) * (1 - gz))
    out.scatter_add_(-1, idx, out1)
    # - corner 101
    idx = sub2ind_list([gx1, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0, signz1])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * (gx * (1 - gy) * gz)
    out.scatter_add_(-1, idx, out1)
    # - corner 110
    idx = sub2ind_list([gx1, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1, signz0])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * (gx * gy * (1 - gz))
    out.scatter_add_(-1, idx, out1)
    # - corner 111
    idx = sub2ind_list([gx1, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1, signz1])
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * (gx * gy * gz)
    out.scatter_add_(-1, idx, out1)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def grad3d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, iZ) tensor
    g: (B, oX, oY, oZ, 3) tensor
    bound: List{3}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY, oZ, 3) tensor
    """
    dim = 3
    boundx, boundy, boundz = bound
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = torch.unbind(g, -1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)
    gz, gz0, gz1, signz0, signz1 = get_weights_and_indices(gz, nz, boundz)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    out = torch.empty([batch, channel] + list(g.shape[-2:]),
                      dtype=inp.dtype, device=inp.device)
    outx, outy, outz = out.unbind(-1)
    # - corner 000
    idx = sub2ind_list([gx0, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    torch.gather(inp, -1, idx, out=outx)
    outy.copy_(outx)
    outz.copy_(outx)
    sign = make_sign([signx0, signy0, signz0])
    if sign is not None:
        out *= sign.unsqueeze(-1)
    outx *= - (1 - gy) * (1 - gz)
    outy *= - (1 - gx) * (1 - gz)
    outz *= - (1 - gx) * (1 - gy)
    # - corner 001
    idx = sub2ind_list([gx0, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy0, signz1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, - (1 - gy) * gz)
    outy.addcmul_(out1, - (1 - gx) * gz)
    outz.addcmul_(out1,   (1 - gx) * (1 - gy))
    # - corner 010
    idx = sub2ind_list([gx0, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz0])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, - gy * (1 - gz))
    outy.addcmul_(out1, (1 - gx) * (1 - gz))
    outz.addcmul_(out1, - (1 - gx) * gy)
    # - corner 011
    idx = sub2ind_list([gx0, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, - gy * gz)
    outy.addcmul_(out1, (1 - gx) * gz)
    outz.addcmul_(out1, (1 - gx) * gy)
    # - corner 100
    idx = sub2ind_list([gx1, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz0])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, (1 - gy) * (1 - gz))
    outy.addcmul_(out1, - gx * (1 - gz))
    outz.addcmul_(out1, - gx * (1 - gy))
    # - corner 101
    idx = sub2ind_list([gx1, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, (1 - gy) * gz)
    outy.addcmul_(out1, - gx * gz)
    outz.addcmul_(out1, gx * (1 - gy))
    # - corner 110
    idx = sub2ind_list([gx1, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz0])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, gy * (1 - gz))
    outy.addcmul_(out1, gx * (1 - gz))
    outz.addcmul_(out1, - gx * gy)
    # - corner 111
    idx = sub2ind_list([gx1, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, gy * gz)
    outy.addcmul_(out1, gx * gz)
    outz.addcmul_(out1, gx * gy)

    if mask is not None:
        out *= mask.unsqueeze(-1)
    out = out.reshape(list(out.shape[:2]) + oshape + [3])
    return out


@torch.jit.script
def pushgrad3d(inp, g, shape: Optional[List[int]], bound: List[Bound],
               extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, iZ, 3) tensor
    g: (B, iX, iY, iZ, 3) tensor
    shape: List{3}[int], optional
    bound: List{3}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 3
    boundx, boundy, boundz = bound
    if inp.shape[-dim-1:-1] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = list(inp.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = g.unbind(-1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1, dim])
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)
    gz, gz0, gz1, signz0, signz1 = get_weights_and_indices(gz, nz, boundz)

    # scatter
    out = torch.zeros([batch, channel, nx*ny*nz],
                      dtype=inp.dtype, device=inp.device)
    # - corner 000
    idx = sub2ind_list([gx0, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0, signz0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= - (1 - gy) * (1 - gz)
    out1y *= - (1 - gx) * (1 - gz)
    out1z *= - (1 - gx) * (1 - gy)
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 001
    idx = sub2ind_list([gx0, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0, signz1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= - (1 - gy) * gz
    out1y *= - (1 - gx) * gz
    out1z *= (1 - gx) * (1 - gy)
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 010
    idx = sub2ind_list([gx0, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1, signz0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= - gy * (1 - gz)
    out1y *= (1 - gx) * (1 - gz)
    out1z *= - (1 - gx) * gy
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 011
    idx = sub2ind_list([gx0, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1, signz1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= - gy * gz
    out1y *= (1 - gx) * gz
    out1z *= (1 - gx) * gy
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 100
    idx = sub2ind_list([gx1, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0, signz0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= (1 - gy) * (1 - gz)
    out1y *= - gx * (1 - gz)
    out1z *= - gx * (1 - gy)
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 101
    idx = sub2ind_list([gx1, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0, signz1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= (1 - gy) * gz
    out1y *= - gx * gz
    out1z *= gx * (1 - gy)
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 110
    idx = sub2ind_list([gx1, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1, signz0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= gy * (1 - gz)
    out1y *= gx * (1 - gz)
    out1z *= - gx * gy
    out.scatter_add_(-1, idx, out1x + out1y + out1z)
    # - corner 111
    idx = sub2ind_list([gx1, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1, signz1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y, out1z = out1.unbind(-1)
    out1x *= gy * gz
    out1y *= gx * gz
    out1z *= gx * gy
    out.scatter_add_(-1, idx, out1x + out1y + out1z)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def hess3d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, iZ) tensor
    g: (B, oX, oY, oZ, 3) tensor
    bound: List{3}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY, oZ, 3, 3) tensor
    """
    dim = 3
    boundx, boundy, boundz = bound
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy, gz = torch.unbind(g, -1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny, nz = shape

    # mask of inbounds voxels
    mask = inbounds_mask_3d(extrapolate, gx, gy, gz, nx, ny, nz)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)
    gz, gz0, gz1, signz0, signz1 = get_weights_and_indices(gz, nz, boundz)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    out = torch.empty([batch, channel, g.shape[-2], dim, dim],
                      dtype=inp.dtype, device=inp.device)
    outx, outy, outz = out.unbind(-1)
    outxx, outyx, outzx = outx.unbind(-1)
    outxy, outyy, outzy = outy.unbind(-1)
    outxz, outyz, outzz = outz.unbind(-1)
    # - corner 000
    idx = sub2ind_list([gx0, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    torch.gather(inp, -1, idx, out=outxy)
    outxz.copy_(outxy)
    outyz.copy_(outxy)
    outxx.zero_()
    outyy.zero_()
    outzz.zero_()
    sign = make_sign([signx0, signy0, signz0])
    if sign is not None:
        out *= sign.unsqueeze(-1).unsqueeze(-1)
    outxy *= (1 - gz)
    outxz *= (1 - gy)
    outyz *= (1 - gx)
    # - corner 001
    idx = sub2ind_list([gx0, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy0, signz1])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, gz)
    outxz.addcmul_(out1, - (1 - gy))
    outyz.addcmul_(out1, - (1 - gx))
    # - corner 010
    idx = sub2ind_list([gx0, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz0])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, - (1 - gz))
    outxz.addcmul_(out1, gy)
    outyz.addcmul_(out1, - (1 - gx))
    # - corner 011
    idx = sub2ind_list([gx0, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1, signz1])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, - gz)
    outxz.addcmul_(out1, - gy)
    outyz.addcmul_(out1, (1 - gx))
    # - corner 100
    idx = sub2ind_list([gx1, gy0, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz0])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, - (1 - gz))
    outxz.addcmul_(out1, - (1 - gy))
    outyz.addcmul_(out1, gx)
    # - corner 101
    idx = sub2ind_list([gx1, gy0, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0, signz1])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, - gz)
    outxz.addcmul_(out1, (1 - gy))
    outyz.addcmul_(out1, - gx)
    # - corner 110
    idx = sub2ind_list([gx1, gy1, gz0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz0])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, (1 - gz))
    outxz.addcmul_(out1, - gy)
    outyz.addcmul_(out1, - gx)
    # - corner 111
    idx = sub2ind_list([gx1, gy1, gz1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1, signz1])
    if sign is not None:
        out1 *= sign
    outxy.addcmul_(out1, gz)
    outxz.addcmul_(out1, gy)
    outyz.addcmul_(out1, gx)

    outyx.copy_(outxy)
    outzx.copy_(outxz)
    outzy.copy_(outyz)

    if mask is not None:
        out *= mask.unsqueeze(-1).unsqueeze(-1)
    out = out.reshape(list(out.shape[:2]) + oshape + [dim, dim])
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
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = g.unbind(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny = shape
    
    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    # - corner 00
    idx = sub2ind_list([gx0, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = make_sign([signx0, signy0])
    if sign is not None:
        out = out * sign
    out = out * ((1 - gx) * (1 - gy))
    # - corner 01
    idx = sub2ind_list([gx0, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * ((1 - gx) * gy)
    out = out + out1
    # - corner 10
    idx = sub2ind_list([gx1, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * (1 - gy))
    out = out + out1
    # - corner 11
    idx = sub2ind_list([gx1, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1])
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * (gx * gy)
    out = out + out1

    if mask is not None:
        out *= mask
    out = out.reshape(list(out.shape[:2]) + oshape)
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
    ishape = list(inp.shape[-dim:])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = torch.unbind(g, -1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)

    # scatter
    out = torch.zeros([batch, channel, nx*ny],
                      dtype=inp.dtype, device=inp.device)
    # - corner 00
    idx = sub2ind_list([gx0, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0])
    if sign is not None:
        out1 *= sign
    if mask is not None:
        out1 *= mask
    out1 *= (1 - gx) * (1 - gy)
    out.scatter_add_(-1, idx, out1)
    # - corner 01
    idx = sub2ind_list([gx0, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1])
    if sign is not None:
        out1 *= sign
    if mask is not None:
        out1 *= mask
    out1 *= (1 - gx) * gy
    out.scatter_add_(-1, idx, out1)
    # - corner 10
    idx = sub2ind_list([gx1, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0])
    if sign is not None:
        out1 *= sign
    if mask is not None:
        out1 *= mask
    out1 *= gx * (1 - gy)
    out.scatter_add_(-1, idx, out1)
    # - corner 11
    idx = sub2ind_list([gx1, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1])
    if sign is not None:
        out1 *= sign
    if mask is not None:
        out1 *= mask
    out1 *= gx * gy
    out.scatter_add_(-1, idx, out1)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def grad2d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY) tensor
    g: (B, oX, oY, 2) tensor
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY, 2) tensor
    """
    dim = 2
    boundx, boundy = bound
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = torch.unbind(g, -1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    out = torch.empty([batch, channel] + list(g.shape[-2:]),
                      dtype=inp.dtype, device=inp.device)
    outx, outy = out.unbind(-1)
    # - corner 00
    idx = sub2ind_list([gx0, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    torch.gather(inp, -1, idx, out=outx)
    outy.copy_(outx)
    sign = make_sign([signx0, signy0])
    if sign is not None:
        out *= sign.unsqueeze(-1)
    outx *= - (1 - gy)
    outy *= - (1 - gx)
    # - corner 01
    idx = sub2ind_list([gx0, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, - gy)
    outy.addcmul_(out1, (1 - gx))
    # - corner 10
    idx = sub2ind_list([gx1, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, (1 - gy))
    outy.addcmul_(out1, - gx)
    # - corner 11
    idx = sub2ind_list([gx1, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1])
    if sign is not None:
        out1 *= sign
    outx.addcmul_(out1, gy)
    outy.addcmul_(out1, gx)

    if mask is not None:
        out *= mask.unsqueeze(-1)
    out = out.reshape(list(out.shape[:2]) + oshape + [dim])
    return out


@torch.jit.script
def pushgrad2d(inp, g, shape: Optional[List[int]], bound: List[Bound],
               extrapolate: int = 1):
    """
    inp: (B, C, iX, iY, 2) tensor
    g: (B, iX, iY, 2) tensor
    shape: List{2}[int], optional
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 2
    boundx, boundy = bound
    if inp.shape[-dim-1:-1] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = list(inp.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = g.unbind(-1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1, dim])
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)

    # scatter
    out = torch.zeros([batch, channel, nx*ny],
                      dtype=inp.dtype, device=inp.device)
    # - corner 00
    idx = sub2ind_list([gx0, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y = out1.unbind(-1)
    out1x *= - (1 - gy)
    out1y *= - (1 - gx)
    out.scatter_add_(-1, idx, out1x + out1y)
    # - corner 01
    idx = sub2ind_list([gx0, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx0, signy1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y = out1.unbind(-1)
    out1x *= - gy
    out1y *= (1 - gx)
    out.scatter_add_(-1, idx, out1x + out1y)
    # - corner 10
    idx = sub2ind_list([gx1, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy0])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y = out1.unbind(-1)
    out1x *= (1 - gy)
    out1y *= - gx
    out.scatter_add_(-1, idx, out1x + out1y)
    # - corner 11
    idx = sub2ind_list([gx1, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = make_sign([signx1, signy1])
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x, out1y = out1.unbind(-1)
    out1x *= gy
    out1y *= gx
    out.scatter_add_(-1, idx, out1x + out1y)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def hess2d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX, iY) tensor
    g: (B, oX, oY, 2) tensor
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, oY, 2, 2) tensor
    """
    dim = 2
    boundx, boundy = bound
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx, gy = torch.unbind(g, -1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx, ny = shape

    # mask of inbounds voxels
    mask = inbounds_mask_2d(extrapolate, gx, gy, nx, ny)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)
    gy, gy0, gy1, signy0, signy1 = get_weights_and_indices(gy, ny, boundy)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    out = torch.empty([batch, channel, g.shape[-2], dim, dim],
                      dtype=inp.dtype, device=inp.device)
    outx, outy = out.unbind(-1)
    outxx, outyx = outx.unbind(-1)
    outxy, outyy = outy.unbind(-1)
    # - corner 00
    idx = sub2ind_list([gx0, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    torch.gather(inp, -1, idx, out=outxy)
    outxx.zero_()
    outyy.zero_()
    sign = make_sign([signx0, signy0])
    if sign is not None:
        out *= sign.unsqueeze(-1).unsqueeze(-1)
    outxy *= 1
    # - corner 01
    idx = sub2ind_list([gx0, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx0, signy1])
    if sign is not None:
        out1 *= sign
    outxy.add_(out1, alpha=-1)
    # - corner 10
    idx = sub2ind_list([gx1, gy0], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy0])
    if sign is not None:
        out1 *= sign
    outxy.add_(out1, alpha=-1)
    # - corner 11
    idx = sub2ind_list([gx1, gy1], shape)
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = make_sign([signx1, signy1])
    if sign is not None:
        out1 *= sign
    outxy.add_(out1)

    outyx.copy_(outxy)

    if mask is not None:
        out *= mask.unsqueeze(-1).unsqueeze(-1)
    out = out.reshape(list(out.shape[:2]) + oshape + [dim, dim])
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
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    # - corner 0
    idx = gx0
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out = inp.gather(-1, idx)
    sign = signx0
    if sign is not None:
        out = out * sign
    out = out * (1 - gx)
    # - corner 1
    idx = gx1
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = signx1
    if sign is not None:
        out1 = out1 * sign
    out1 = out1 * gx
    out = out + out1

    if mask is not None:
        out *= mask
    out = out.reshape(list(out.shape[:2]) + oshape)
    return out


@torch.jit.script
def push1d(inp, g, shape: Optional[List[int]], bound: List[Bound],
           extrapolate: int = 1):
    """
    inp: (B, C, iX, iY) tensor
    g: (B, iX, iY, 2) tensor
    shape: List{2}[int], optional
    bound: List{2}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 1
    boundx = bound[0]
    if inp.shape[-dim:] != g.shape[-dim-1:-1]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = list(inp.shape[-dim:])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)

    # scatter
    out = torch.zeros([batch, channel, nx],
                      dtype=inp.dtype, device=inp.device)
    # - corner 0
    idx = gx0
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = signx0
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * (1 - gx)
    out.scatter_add_(-1, idx, out1)
    # - corner 1
    idx = gx1
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = signx1
    if sign is not None:
        out1 = out1 * sign
    if mask is not None:
        out1 = out1 * mask
    out1 = out1 * gx
    out.scatter_add_(-1, idx, out1)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def grad1d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX) tensor
    g: (B, oX, 1) tensor
    bound: List{1}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, 1) tensor
    """
    dim = 1
    boundx = bound[0]
    oshape = list(g.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    batch = max(inp.shape[0], gx.shape[0])
    channel = inp.shape[1]
    shape = list(inp.shape[-dim:])
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)

    # gather
    inp = inp.reshape(list(inp.shape[:2]) + [-1])
    out = torch.empty([batch, channel] + list(g.shape[-2:]),
                      dtype=inp.dtype, device=inp.device)
    outx = out.squeeze(-1)
    # - corner 0
    idx = gx0
    idx = idx.expand([batch, channel, idx.shape[-1]])
    torch.gather(inp, -1, idx, out=outx)
    sign = signx0
    if sign is not None:
        out *= sign.unsqueeze(-1)
    outx.neg_()
    # - corner 1
    idx = gx1
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.gather(-1, idx)
    sign = signx1
    if sign is not None:
        out1 *= sign
    outx.add_(out1)

    if mask is not None:
        out *= mask.unsqueeze(-1)
    out = out.reshape(list(out.shape[:2]) + oshape + [dim])
    return out


@torch.jit.script
def pushgrad1d(inp, g, shape: Optional[List[int]], bound: List[Bound],
               extrapolate: int = 1):
    """
    inp: (B, C, iX, 1) tensor
    g: (B, iX, 1) tensor
    shape: List{1}[int], optional
    bound: List{1}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, *shape) tensor
    """
    dim = 1
    boundx = bound[0]
    if inp.shape[-2] != g.shape[-2]:
        raise ValueError('Input and grid should have the same spatial shape')
    ishape = list(inp.shape[-dim-1:-1])
    g = g.reshape([g.shape[0], 1, -1, dim])
    gx = g.squeeze(-1)
    inp = inp.reshape(list(inp.shape[:2]) + [-1, dim])
    batch = max(inp.shape[0], g.shape[0])
    channel = inp.shape[1]

    if shape is None:
        shape = ishape
    shape = list(shape)
    nx = shape[0]

    # mask of inbounds voxels
    mask = inbounds_mask_1d(extrapolate, gx, nx)

    # corners
    # (upper weight, lower corner, upper corner, lower sign, upper sign)
    gx, gx0, gx1, signx0, signx1 = get_weights_and_indices(gx, nx, boundx)

    # scatter
    out = torch.zeros([batch, channel, nx], dtype=inp.dtype, device=inp.device)
    # - corner 000
    idx = gx0
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = signx0
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x = out1.squeeze(-1)
    out1x.neg_()
    out.scatter_add_(-1, idx, out1x)
    # - corner 100
    idx = gx1
    idx = idx.expand([batch, channel, idx.shape[-1]])
    out1 = inp.clone()
    sign = signx1
    if sign is not None:
        out1 *= sign.unsqueeze(-1)
    if mask is not None:
        out1 *= mask.unsqueeze(-1)
    out1x = out1.squeeze(-1)
    out.scatter_add_(-1, idx, out1x)

    out = out.reshape(list(out.shape[:2]) + shape)
    return out


@torch.jit.script
def hess1d(inp, g, bound: List[Bound], extrapolate: int = 1):
    """
    inp: (B, C, iX) tensor
    g: (B, oX, 1) tensor
    bound: List{1}[Bound] tensor
    extrapolate: ExtrapolateType
    returns: (B, C, oX, 1, 1) tensor
    """
    batch = max(inp.shape[0], g.shape[0])
    return torch.zeros([batch, inp.shape[1], g.shape[1], 1, 1],
                       dtype=inp.dtype, device=inp.device)