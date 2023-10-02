try:
    import jitfields
    available = True
except (ImportError, ModuleNotFoundError):
    jitfields = None
    available = False
from .utils import make_list
import torch


def first2last(input, ndim):
    insert = input.dim() <= ndim
    if insert:
        input = input.unsqueeze(-1)
    else:
        input = torch.movedim(input, -ndim-1, -1)
    return input, insert


def last2first(input, ndim, inserted, grad=False):
    if inserted:
        input = input.squeeze(-1 - grad)
    else:
        input = torch.movedim(input, -1 - grad, -ndim-1 - grad)
    return input


def grid_pull(input, grid, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    ndim = grid.shape[-1]
    input, inserted = first2last(input, ndim)
    input = jitfields.pull(input, grid, order=interpolation, bound=bound,
                           extrapolate=extrapolate, prefilter=prefilter)
    input = last2first(input, ndim, inserted)
    return input


def grid_push(input, grid, shape=None, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    ndim = grid.shape[-1]
    input, inserted = first2last(input, ndim)
    input = jitfields.push(input, grid, shape, order=interpolation, bound=bound,
                           extrapolate=extrapolate, prefilter=prefilter)
    input = last2first(input, ndim, inserted)
    return input


def grid_count(grid, shape=None, interpolation='linear', bound='zero',
               extrapolate=False):
    return jitfields.count(grid, shape, order=interpolation, bound=bound,
                           extrapolate=extrapolate)


def grid_grad(input, grid, interpolation='linear', bound='zero',
              extrapolate=False, prefilter=False):
    ndim = grid.shape[-1]
    input, inserted = first2last(input, ndim)
    input = jitfields.grad(input, grid, order=interpolation, bound=bound,
                           extrapolate=extrapolate, prefilter=prefilter)
    input = last2first(input, ndim, inserted, True)
    return input


def spline_coeff(input, interpolation='linear', bound='dct2', dim=-1,
                 inplace=False):
    func = jitfields.spline_coeff_ if inplace else jitfields.spline_coeff
    return func(input, interpolation, bound=bound, dim=dim)


def spline_coeff_nd(input, interpolation='linear', bound='dct2', dim=None,
                    inplace=False):
    func = jitfields.spline_coeff_nd_ if inplace else jitfields.spline_coeff_nd
    return func(input, interpolation, bound=bound, ndim=dim)


def resize(image, factor=None, shape=None, anchor='c',
           interpolation=1, prefilter=True, **kwargs):
    kwargs.setdefault('bound', 'nearest')
    ndim = max(len(make_list(factor or [])),
               len(make_list(shape or [])),
               len(make_list(anchor or []))) or (image.dim() - 2)
    return jitfields.resize(image, factor=factor, shape=shape, ndim=ndim,
                            anchor=anchor, order=interpolation,
                            bound=kwargs['bound'], prefilter=prefilter)


def restrict(image, factor=None, shape=None, anchor='c',
             interpolation=1, reduce_sum=False, **kwargs):
    kwargs.setdefault('bound', 'nearest')
    ndim = max(len(make_list(factor or [])),
               len(make_list(shape or [])),
               len(make_list(anchor or []))) or (image.dim() - 2)
    return jitfields.restrict(image, factor=factor, shape=shape, ndim=ndim,
                              anchor=anchor, order=interpolation,
                              bound=kwargs['bound'], reduce_sum=reduce_sum)
