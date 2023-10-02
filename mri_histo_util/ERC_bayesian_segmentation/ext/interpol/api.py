"""High level interpolation API"""

__all__ = ['grid_pull', 'grid_push', 'grid_count', 'grid_grad',
           'spline_coeff', 'spline_coeff_nd',
           'identity_grid', 'add_identity_grid', 'add_identity_grid_']

import torch
from .utils import expanded_shape, matvec
from .jit_utils import movedim1, meshgrid
from .autograd import (GridPull, GridPush, GridCount, GridGrad,
                       SplineCoeff, SplineCoeffND)
from . import backend, jitfields

_doc_interpolation = \
"""`interpolation` can be an int, a string or an InterpolationType.
    Possible values are:
        - 0 or 'nearest'
        - 1 or 'linear'
        - 2 or 'quadratic'
        - 3 or 'cubic'
        - 4 or 'fourth'
        - 5 or 'fifth'
        - etc.
    A list of values can be provided, in the order [W, H, D],
    to specify dimension-specific interpolation orders."""

_doc_bound = \
"""`bound` can be an int, a string or a BoundType.
    Possible values are:
        - 'replicate'  or 'nearest'     :  a  a  a  |  a  b  c  d  |  d  d  d
        - 'dct1'       or 'mirror'      :  d  c  b  |  a  b  c  d  |  c  b  a
        - 'dct2'       or 'reflect'     :  c  b  a  |  a  b  c  d  |  d  c  b
        - 'dst1'       or 'antimirror'  : -b -a  0  |  a  b  c  d  |  0 -d -c
        - 'dst2'       or 'antireflect' : -c -b -a  |  a  b  c  d  | -d -c -b
        - 'dft'        or 'wrap'        :  b  c  d  |  a  b  c  d  |  a  b  c
        - 'zero'       or 'zeros'       :  0  0  0  |  a  b  c  d  |  0  0  0
    A list of values can be provided, in the order [W, H, D],
    to specify dimension-specific boundary conditions.
    Note that
    - `dft` corresponds to circular padding
    - `dct2` corresponds to Neumann boundary conditions (symmetric)
    - `dst2` corresponds to Dirichlet boundary conditions (antisymmetric)
    See https://en.wikipedia.org/wiki/Discrete_cosine_transform
        https://en.wikipedia.org/wiki/Discrete_sine_transform"""

_doc_bound_coeff = \
"""`bound` can be an int, a string or a BoundType. 
    Possible values are:
        - 'replicate'  or 'nearest'     :  a  a  a  |  a  b  c  d  |  d  d  d
        - 'dct1'       or 'mirror'      :  d  c  b  |  a  b  c  d  |  c  b  a
        - 'dct2'       or 'reflect'     :  c  b  a  |  a  b  c  d  |  d  c  b
        - 'dst1'       or 'antimirror'  : -b -a  0  |  a  b  c  d  |  0 -d -c
        - 'dst2'       or 'antireflect' : -c -b -a  |  a  b  c  d  | -d -c -b
        - 'dft'        or 'wrap'        :  b  c  d  |  a  b  c  d  |  a  b  c
        - 'zero'       or 'zeros'       :  0  0  0  |  a  b  c  d  |  0  0  0
    A list of values can be provided, in the order [W, H, D],
    to specify dimension-specific boundary conditions.
    Note that
    - `dft` corresponds to circular padding
    - `dct1` corresponds to mirroring about the center of the first/last voxel
    - `dct2` corresponds to mirroring about the edge of the first/last voxel
    See https://en.wikipedia.org/wiki/Discrete_cosine_transform
        https://en.wikipedia.org/wiki/Discrete_sine_transform
        
    /!\ Only 'dct1', 'dct2' and 'dft' are implemented for interpolation
        orders >= 6."""

_ref_coeff = \
"""..[1]  M. Unser, A. Aldroubi and M. Eden.
       "B-Spline Signal Processing: Part I-Theory,"
       IEEE Transactions on Signal Processing 41(2):821-832 (1993).
..[2]  M. Unser, A. Aldroubi and M. Eden.
       "B-Spline Signal Processing: Part II-Efficient Design and Applications,"
       IEEE Transactions on Signal Processing 41(2):834-848 (1993).
..[3]  M. Unser.
       "Splines: A Perfect Fit for Signal and Image Processing,"
       IEEE Signal Processing Magazine 16(6):22-38 (1999).
"""


def _preproc(grid, input=None, mode=None):
    """Preprocess tensors for pull/push/count/grad

    Low level bindings expect inputs of shape
    [batch, channel, *spatial] and [batch, *spatial, dim], whereas
    the high level python API accepts inputs of shape
    [..., [channel], *spatial] and [..., *spatial, dim].

    This function broadcasts and reshapes the input tensors accordingly.
            /!\\ This *can* trigger large allocations /!\\
    """
    dim = grid.shape[-1]
    if input is None:
        spatial = grid.shape[-dim-1:-1]
        batch = grid.shape[:-dim-1]
        grid = grid.reshape([-1, *spatial, dim])
        info = dict(batch=batch, channel=[1] if batch else [], dim=dim)
        return grid, info

    grid_spatial = grid.shape[-dim-1:-1]
    grid_batch = grid.shape[:-dim-1]
    input_spatial = input.shape[-dim:]
    channel = 0 if input.dim() == dim else input.shape[-dim-1]
    input_batch = input.shape[:-dim-1]

    if mode == 'push':
        grid_spatial = input_spatial = expanded_shape(grid_spatial, input_spatial)

    # broadcast and reshape
    batch = expanded_shape(grid_batch, input_batch)
    grid = grid.expand([*batch, *grid_spatial, dim])
    grid = grid.reshape([-1, *grid_spatial, dim])
    input = input.expand([*batch, channel or 1, *input_spatial])
    input = input.reshape([-1, channel or 1, *input_spatial])

    out_channel = [channel] if channel else ([1] if batch else [])
    info = dict(batch=batch, channel=out_channel, dim=dim)
    return grid, input, info


def _postproc(out, shape_info, mode):
    """Postprocess tensors for pull/push/count/grad"""
    dim = shape_info['dim']
    if mode != 'grad':
        spatial = out.shape[-dim:]
        feat = []
    else:
        spatial = out.shape[-dim-1:-1]
        feat = [out.shape[-1]]
    batch = shape_info['batch']
    channel = shape_info['channel']

    out = out.reshape([*batch, *channel, *spatial, *feat])
    return out


def grid_pull(input, grid, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    """Sample an image with respect to a deformation field.

    Notes
    -----
    {interpolation}

    {bound}
    
    If the input dtype is not a floating point type, the input image is 
    assumed to contain labels. Then, unique labels are extracted 
    and resampled individually, making them soft labels. Finally, 
    the label map is reconstructed from the individual soft labels by 
    assigning the label with maximum soft value.

    Parameters
    ----------
    input : (..., [channel], *inshape) tensor
        Input image.
    grid : (..., *outshape, dim) tensor
        Transformation field.
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType or sequence[BoundType], default='zero'
        Boundary conditions.
    extrapolate : bool or int, default=True
        Extrapolate out-of-bound data.
    prefilter : bool, default=False
        Apply spline pre-filter (= interpolates the input)

    Returns
    -------
    output : (..., [channel], *outshape) tensor
        Deformed image.

    """
    if backend.jitfields and jitfields.available:
        return jitfields.grid_pull(input, grid, interpolation, bound,
                                   extrapolate, prefilter)

    grid, input, shape_info = _preproc(grid, input)
    batch, channel = input.shape[:2]
    dim = grid.shape[-1]

    if not input.dtype.is_floating_point:
        # label map -> specific processing
        out = input.new_zeros([batch, channel, *grid.shape[1:-1]])
        pmax = grid.new_zeros([batch, channel, *grid.shape[1:-1]])
        for label in input.unique():
            soft = (input == label).to(grid.dtype)
            if prefilter:
                input = spline_coeff_nd(soft, interpolation=interpolation,
                                        bound=bound, dim=dim, inplace=True)
            soft = GridPull.apply(soft, grid, interpolation, bound, extrapolate)
            out[soft > pmax] = label
            pmax = torch.max(pmax, soft)
    else:
        if prefilter:
            input = spline_coeff_nd(input, interpolation=interpolation,
                                    bound=bound, dim=dim)
        out = GridPull.apply(input, grid, interpolation, bound, extrapolate)

    return _postproc(out, shape_info, mode='pull')


def grid_push(input, grid, shape=None, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    """Splat an image with respect to a deformation field (pull adjoint).

    Notes
    -----
    {interpolation}

    {bound}

    Parameters
    ----------
    input : (..., [channel], *inshape) tensor
        Input image.
    grid : (..., *inshape, dim) tensor
        Transformation field.
    shape : sequence[int], default=inshape
        Output shape
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType, or sequence[BoundType], default='zero'
        Boundary conditions.
    extrapolate : bool or int, default=True
        Extrapolate out-of-bound data.
    prefilter : bool, default=False
        Apply spline pre-filter.

    Returns
    -------
    output : (..., [channel], *shape) tensor
        Spatted image.

    """
    if backend.jitfields and jitfields.available:
        return jitfields.grid_push(input, grid, shape, interpolation, bound,
                                   extrapolate, prefilter)

    grid, input, shape_info = _preproc(grid, input, mode='push')
    dim = grid.shape[-1]

    if shape is None:
        shape = tuple(input.shape[2:])

    out = GridPush.apply(input, grid, shape, interpolation, bound, extrapolate)
    if prefilter:
        out = spline_coeff_nd(out, interpolation=interpolation, bound=bound,
                              dim=dim, inplace=True)
    return _postproc(out, shape_info, mode='push')


def grid_count(grid, shape=None, interpolation='linear', bound='zero',
               extrapolate=False):
    """Splatting weights with respect to a deformation field (pull adjoint).

    Notes
    -----
    {interpolation}

    {bound}

    Parameters
    ----------
    grid : (..., *inshape, dim) tensor
        Transformation field.
    shape : sequence[int], default=inshape
        Output shape
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType, or sequence[BoundType], default='zero'
        Boundary conditions.
    extrapolate : bool or int, default=True
        Extrapolate out-of-bound data.

    Returns
    -------
    output : (..., [1], *shape) tensor
        Splatted weights.

    """
    if backend.jitfields and jitfields.available:
        return jitfields.grid_count(grid, shape, interpolation, bound, extrapolate)

    grid, shape_info = _preproc(grid)
    out = GridCount.apply(grid, shape, interpolation, bound, extrapolate)
    return _postproc(out, shape_info, mode='count')


def grid_grad(input, grid, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    """Sample spatial gradients of an image with respect to a deformation field.
    
    Notes
    -----
    {interpolation}

    {bound}

    Parameters
    ----------
    input : (..., [channel], *inshape) tensor
        Input image.
    grid : (..., *inshape, dim) tensor
        Transformation field.
    shape : sequence[int], default=inshape
        Output shape
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType, or sequence[BoundType], default='zero'
        Boundary conditions.
    extrapolate : bool or int, default=True
        Extrapolate out-of-bound data.
    prefilter : bool, default=False
        Apply spline pre-filter (= interpolates the input)

    Returns
    -------
    output : (..., [channel], *shape, dim) tensor
        Sampled gradients.

    """
    if backend.jitfields and jitfields.available:
        return jitfields.grid_grad(input, grid, interpolation, bound,
                                   extrapolate, prefilter)

    grid, input, shape_info = _preproc(grid, input)
    dim = grid.shape[-1]
    if prefilter:
        input = spline_coeff_nd(input, interpolation, bound, dim)
    out = GridGrad.apply(input, grid, interpolation, bound, extrapolate)
    return _postproc(out, shape_info, mode='grad')


def spline_coeff(input, interpolation='linear', bound='dct2', dim=-1,
                 inplace=False):
    """Compute the interpolating spline coefficients, for a given spline order
    and boundary conditions, along a single dimension.

    Notes
    -----
    {interpolation}

    {bound}

    References
    ----------
    {ref}


    Parameters
    ----------
    input : tensor
        Input image.
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType or sequence[BoundType], default='dct1'
        Boundary conditions.
    dim : int, default=-1
        Dimension along which to process
    inplace : bool, default=False
        Process the volume in place.

    Returns
    -------
    output : tensor
        Coefficient image.

    """
    # This implementation is based on the file bsplines.c in SPM12, written
    # by John Ashburner, which is itself based on the file coeff.c,
    # written by Philippe Thevenaz: http://bigwww.epfl.ch/thevenaz/interpolation
    # . DCT1 boundary conditions were derived by Thevenaz and Unser.
    # . DFT boundary conditions were derived by John Ashburner.
    # SPM12 is released under the GNU-GPL v2 license.
    # Philippe Thevenaz's code does not have an explicit license as far
    # as we know.
    if backend.jitfields and jitfields.available:
        return jitfields.spline_coeff(input, interpolation, bound,
                                      dim, inplace)

    out = SplineCoeff.apply(input, bound, interpolation, dim, inplace)
    return out


def spline_coeff_nd(input, interpolation='linear', bound='dct2', dim=None,
                    inplace=False):
    """Compute the interpolating spline coefficients, for a given spline order
    and boundary conditions, along the last `dim` dimensions.

    Notes
    -----
    {interpolation}

    {bound}

    References
    ----------
    {ref}

    Parameters
    ----------
    input : (..., *spatial) tensor
        Input image.
    interpolation : int or sequence[int], default=1
        Interpolation order.
    bound : BoundType or sequence[BoundType], default='dct1'
        Boundary conditions.
    dim : int, default=-1
        Number of spatial dimensions
    inplace : bool, default=False
        Process the volume in place.

    Returns
    -------
    output : (..., *spatial) tensor
        Coefficient image.

    """
    # This implementation is based on the file bsplines.c in SPM12, written
    # by John Ashburner, which is itself based on the file coeff.c,
    # written by Philippe Thevenaz: http://bigwww.epfl.ch/thevenaz/interpolation
    # . DCT1 boundary conditions were derived by Thevenaz and Unser.
    # . DFT boundary conditions were derived by John Ashburner.
    # SPM12 is released under the GNU-GPL v2 license.
    # Philippe Thevenaz's code does not have an explicit license as far
    # as we know.
    if backend.jitfields and jitfields.available:
        return jitfields.spline_coeff_nd(input, interpolation, bound,
                                         dim, inplace)

    out = SplineCoeffND.apply(input, bound, interpolation, dim, inplace)
    return out


grid_pull.__doc__ = grid_pull.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound)
grid_push.__doc__ = grid_push.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound)
grid_count.__doc__ = grid_count.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound)
grid_grad.__doc__ = grid_grad.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound)
spline_coeff.__doc__ = spline_coeff.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound_coeff, ref=_ref_coeff)
spline_coeff_nd.__doc__ = spline_coeff_nd.__doc__.format(
    interpolation=_doc_interpolation, bound=_doc_bound_coeff, ref=_ref_coeff)

# aliases
pull = grid_pull
push = grid_push
count = grid_count


def identity_grid(shape, dtype=None, device=None):
    """Returns an identity deformation field.

    Parameters
    ----------
    shape : (dim,) sequence of int
        Spatial dimension of the field.
    dtype : torch.dtype, default=`get_default_dtype()`
        Data type.
    device torch.device, optional
        Device.

    Returns
    -------
    grid : (*shape, dim) tensor
        Transformation field

    """
    mesh1d = [torch.arange(float(s), dtype=dtype, device=device)
              for s in shape]
    grid = torch.stack(meshgrid(mesh1d), dim=-1)
    return grid


@torch.jit.script
def add_identity_grid_(disp):
    """Adds the identity grid to a displacement field, inplace.

    Parameters
    ----------
    disp : (..., *spatial, dim) tensor
        Displacement field

    Returns
    -------
    grid : (..., *spatial, dim) tensor
        Transformation field

    """
    dim = disp.shape[-1]
    spatial = disp.shape[-dim-1:-1]
    mesh1d = [torch.arange(s, dtype=disp.dtype, device=disp.device)
              for s in spatial]
    grid = meshgrid(mesh1d)
    disp = movedim1(disp, -1, 0)
    for i, grid1 in enumerate(grid):
        disp[i].add_(grid1)
    disp = movedim1(disp, 0, -1)
    return disp


@torch.jit.script
def add_identity_grid(disp):
    """Adds the identity grid to a displacement field.

    Parameters
    ----------
    disp : (..., *spatial, dim) tensor
        Displacement field

    Returns
    -------
    grid : (..., *spatial, dim) tensor
        Transformation field

    """
    return add_identity_grid_(disp.clone())


def affine_grid(mat, shape):
    """Create a dense transformation grid from an affine matrix.

    Parameters
    ----------
    mat : (..., D[+1], D+1) tensor
        Affine matrix (or matrices).
    shape : (D,) sequence[int]
        Shape of the grid, with length D.

    Returns
    -------
    grid : (..., *shape, D) tensor
        Dense transformation grid

    """
    mat = torch.as_tensor(mat)
    shape = list(shape)
    nb_dim = mat.shape[-1] - 1
    if nb_dim != len(shape):
        raise ValueError('Dimension of the affine matrix ({}) and shape ({}) '
                         'are not the same.'.format(nb_dim, len(shape)))
    if mat.shape[-2] not in (nb_dim, nb_dim+1):
        raise ValueError('First argument should be matrces of shape '
                         '(..., {0}, {1}) or (..., {1], {1}) but got {2}.'
                         .format(nb_dim, nb_dim+1, mat.shape))
    batch_shape = mat.shape[:-2]
    grid = identity_grid(shape, mat.dtype, mat.device)
    if batch_shape:
        for _ in range(len(batch_shape)):
            grid = grid.unsqueeze(0)
        for _ in range(nb_dim):
            mat = mat.unsqueeze(-1)
    lin = mat[..., :nb_dim, :nb_dim]
    off = mat[..., :nb_dim, -1]
    grid = matvec(lin, grid) + off
    return grid
