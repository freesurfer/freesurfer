__all__ = ['restrict']

from .api import grid_push
from .utils import make_list, meshgrid_ij
from . import backend, jitfields
import torch


def restrict(image, factor=None, shape=None, anchor='c',
             interpolation=1, reduce_sum=False, **kwargs):
    """Restrict an image by a factor or to a specific shape.

    Notes
    -----
    .. A least one of `factor` and `shape` must be specified
    .. If `anchor in ('centers', 'edges')`, exactly one of `factor` or
       `shape must be specified.
    .. If `anchor in ('first', 'last')`, `factor` must be provided even
       if `shape` is specified.
    .. Because of rounding, it is in general not assured that
       `resize(resize(x, f), 1/f)` returns a tensor with the same shape as x.

        edges          centers          first           last
    e - + - + - e   + - + - + - +   + - + - + - +   + - + - + - +
    | . | . | . |   | c | . | c |   | f | . | . |   | . | . | . |
    + _ + _ + _ +   + _ + _ + _ +   + _ + _ + _ +   + _ + _ + _ +
    | . | . | . |   | . | . | . |   | . | . | . |   | . | . | . |
    + _ + _ + _ +   + _ + _ + _ +   + _ + _ + _ +   + _ + _ + _ +
    | . | . | . |   | c | . | c |   | . | . | . |   | . | . | l |
    e _ + _ + _ e   + _ + _ + _ +   + _ + _ + _ +   + _ + _ + _ +

    Parameters
    ----------
    image : (batch, channel, *inshape) tensor
        Image to resize
    factor : float or list[float], optional
        Resizing factor
        * > 1 : larger image <-> smaller voxels
        * < 1 : smaller image <-> larger voxels
    shape : (ndim,) list[int], optional
        Output shape
    anchor : {'centers', 'edges', 'first', 'last'} or list, default='centers'
        * In cases 'c' and 'e', the volume shape is multiplied by the
          zoom factor (and eventually truncated), and two anchor points
          are used to determine the voxel size.
        * In cases 'f' and 'l', a single anchor point is used so that
          the voxel size is exactly divided by the zoom factor.
          This case with an integer factor corresponds to subslicing
          the volume (e.g., `vol[::f, ::f, ::f]`).
        * A list of anchors (one per dimension) can also be provided.
    interpolation : int or sequence[int], default=1
        Interpolation order.
    reduce_sum : bool, default=False
        Do not normalize by the number of accumulated values per voxel

    Returns
    -------
    restricted : (batch, channel, *shape) tensor
        Restricted image

    """
    if backend.jitfields and jitfields.available:
        return jitfields.restrict(image, factor, shape, anchor,
                                  interpolation, reduce_sum, **kwargs)

    factor = make_list(factor) if factor else []
    shape = make_list(shape) if shape else []
    anchor = make_list(anchor)
    nb_dim = max(len(factor), len(shape), len(anchor)) or (image.dim() - 2)
    anchor = [a[0].lower() for a in make_list(anchor, nb_dim)]
    bck = dict(dtype=image.dtype, device=image.device)

    # compute output shape
    inshape = image.shape[-nb_dim:]
    if factor:
        factor = make_list(factor, nb_dim)
    elif not shape:
        raise ValueError('One of `factor` or `shape` must be provided')
    if shape:
        shape = make_list(shape, nb_dim)
    else:
        shape = [int(i/f) for i, f in zip(inshape, factor)]

    if not factor:
        factor = [i/o for o, i in zip(shape, inshape)]

    # compute transformation grid
    lin = []
    fullscale = 1
    for anch, f, inshp, outshp in zip(anchor, factor, inshape, shape):
        if anch == 'c':    # centers
            lin.append(torch.linspace(0, outshp - 1, inshp, **bck))
            fullscale *= (inshp - 1) / (outshp - 1)
        elif anch == 'e':  # edges
            scale = outshp / inshp
            shift = 0.5 * (scale - 1)
            fullscale *= scale
            lin.append(torch.arange(0., inshp, **bck) * scale + shift)
        elif anch == 'f':  # first voxel
            # scale = 1/f
            # shift = 0
            fullscale *= 1/f
            lin.append(torch.arange(0., inshp, **bck) / f)
        elif anch == 'l':  # last voxel
            # scale = 1/f
            shift = (outshp - 1) - (inshp - 1) / f
            fullscale *= 1/f
            lin.append(torch.arange(0., inshp, **bck) / f + shift)
        else:
            raise ValueError('Unknown anchor {}'.format(anch))

    # scatter
    kwargs.setdefault('bound', 'nearest')
    kwargs.setdefault('extrapolate', True)
    kwargs.setdefault('interpolation', interpolation)
    kwargs.setdefault('prefilter', False)
    grid = torch.stack(meshgrid_ij(*lin), dim=-1)
    resized = grid_push(image, grid, shape, **kwargs)
    if not reduce_sum:
        resized /= fullscale

    return resized
