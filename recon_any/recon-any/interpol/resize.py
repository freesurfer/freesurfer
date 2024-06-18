"""
Resize functions (equivalent to scipy's zoom, pytorch's interpolate)
based on grid_pull.
"""
__all__ = ['resize']

from .api import grid_pull
from .utils import make_list, meshgrid_ij
from . import backend, jitfields
import torch


def resize(image, factor=None, shape=None, anchor='c',
           interpolation=1, prefilter=True, **kwargs):
    """Resize an image by a factor or to a specific shape.

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
    prefilter : bool, default=True
        Apply spline pre-filter (= interpolates the input)

    Returns
    -------
    resized : (batch, channel, *shape) tensor
        Resized image

    """
    if backend.jitfields and jitfields.available:
        return jitfields.resize(image, factor, shape, anchor,
                                interpolation, prefilter, **kwargs)

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
        shape = [int(i*f) for i, f in zip(inshape, factor)]

    if not factor:
        factor = [o/i for o, i in zip(shape, inshape)]

    # compute transformation grid
    lin = []
    for anch, f, inshp, outshp in zip(anchor, factor, inshape, shape):
        if anch == 'c':    # centers
            lin.append(torch.linspace(0, inshp - 1, outshp, **bck))
        elif anch == 'e':  # edges
            scale = inshp / outshp
            shift = 0.5 * (scale - 1)
            lin.append(torch.arange(0., outshp, **bck) * scale + shift)
        elif anch == 'f':  # first voxel
            # scale = 1/f
            # shift = 0
            lin.append(torch.arange(0., outshp, **bck) / f)
        elif anch == 'l':  # last voxel
            # scale = 1/f
            shift = (inshp - 1) - (outshp - 1) / f
            lin.append(torch.arange(0., outshp, **bck) / f + shift)
        else:
            raise ValueError('Unknown anchor {}'.format(anch))

    # interpolate
    kwargs.setdefault('bound', 'nearest')
    kwargs.setdefault('extrapolate', True)
    kwargs.setdefault('interpolation', interpolation)
    kwargs.setdefault('prefilter', prefilter)
    grid = torch.stack(meshgrid_ij(*lin), dim=-1)
    resized = grid_pull(image, grid, **kwargs)

    return resized

