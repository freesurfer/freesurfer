"""Compute spline interpolating coefficients

These functions are ported from the C routines in SPM's bsplines.c
by John Ashburner, which are themselves ports from Philippe Thevenaz's
code. JA furthermore derived the initial conditions for the DFT ("wrap around")
boundary conditions.

Note that similar routines are available in scipy with boundary conditions
DCT1 ("mirror"), DCT2 ("reflect") and DFT ("wrap"); all derived by P. Thevenaz,
according to the comments. Our DCT2 boundary conditions are ported from
scipy.

Only boundary conditions DCT1, DCT2 and DFT are implemented.

References
----------
..[1]  M. Unser, A. Aldroubi and M. Eden.
       "B-Spline Signal Processing: Part I-Theory,"
       IEEE Transactions on Signal Processing 41(2):821-832 (1993).
..[2]  M. Unser, A. Aldroubi and M. Eden.
       "B-Spline Signal Processing: Part II-Efficient Design and Applications,"
       IEEE Transactions on Signal Processing 41(2):834-848 (1993).
..[3]  M. Unser.
       "Splines: A Perfect Fit for Signal and Image Processing,"
       IEEE Signal Processing Magazine 16(6):22-38 (1999).
"""
import torch
import math
from typing import List, Optional
from .jit_utils import movedim1
from .pushpull import pad_list_int


@torch.jit.script
def get_poles(order: int) -> List[float]:
    empty: List[float] = []
    if order in (0, 1):
        return empty
    if order == 2:
        return [math.sqrt(8.) - 3.]
    if order == 3:
        return [math.sqrt(3.) - 2.]
    if order == 4:
        return [math.sqrt(664. - math.sqrt(438976.)) + math.sqrt(304.) - 19.,
                math.sqrt(664. + math.sqrt(438976.)) - math.sqrt(304.) - 19.]
    if order == 5:
        return [math.sqrt(67.5 - math.sqrt(4436.25)) + math.sqrt(26.25) - 6.5,
                math.sqrt(67.5 + math.sqrt(4436.25)) - math.sqrt(26.25) - 6.5]
    if order == 6:
        return [-0.488294589303044755130118038883789062112279161239377608394,
                -0.081679271076237512597937765737059080653379610398148178525368,
                -0.00141415180832581775108724397655859252786416905534669851652709]
    if order == 7:
        return [-0.5352804307964381655424037816816460718339231523426924148812,
                -0.122554615192326690515272264359357343605486549427295558490763,
                -0.0091486948096082769285930216516478534156925639545994482648003]
    raise NotImplementedError


@torch.jit.script
def get_gain(poles: List[float]) -> float:
    lam: float = 1.
    for pole in poles:
        lam *= (1. - pole) * (1. - 1./pole)
    return lam


@torch.jit.script
def dft_initial(inp, pole: float, dim: int = -1, keepdim: bool = False):

    assert inp.shape[dim] > 1
    max_iter: int = int(math.ceil(-30./math.log(abs(pole))))
    max_iter = min(max_iter, inp.shape[dim])

    poles = torch.as_tensor(pole, dtype=inp.dtype, device=inp.device)
    poles = poles.pow(torch.arange(1, max_iter, dtype=inp.dtype, device=inp.device))
    poles = poles.flip(0)

    inp = movedim1(inp, dim, 0)
    inp0 = inp[0]
    inp = inp[1-max_iter:]
    inp = movedim1(inp, 0, -1)
    out = torch.matmul(inp.unsqueeze(-2), poles.unsqueeze(-1)).squeeze(-1)
    out = out + inp0.unsqueeze(-1)
    if keepdim:
        out = movedim1(out, -1, dim)
    else:
        out = out.squeeze(-1)

    pole = pole ** max_iter
    out = out / (1 - pole)
    return out


@torch.jit.script
def dct1_initial(inp, pole: float, dim: int = -1, keepdim: bool = False):

    n = inp.shape[dim]
    max_iter: int = int(math.ceil(-30./math.log(abs(pole))))

    if max_iter < n:

        poles = torch.as_tensor(pole, dtype=inp.dtype, device=inp.device)
        poles = poles.pow(torch.arange(1, max_iter, dtype=inp.dtype, device=inp.device))

        inp = movedim1(inp, dim, 0)
        inp0 = inp[0]
        inp = inp[1:max_iter]
        inp = movedim1(inp, 0, -1)
        out = torch.matmul(inp.unsqueeze(-2), poles.unsqueeze(-1)).squeeze(-1)
        out = out + inp0.unsqueeze(-1)
        if keepdim:
            out = movedim1(out, -1, dim)
        else:
            out = out.squeeze(-1)

    else:
        max_iter = n

        polen = pole ** (n - 1)
        inp0 = inp[0] + polen * inp[-1]
        inp = inp[1:-1]
        inp = movedim1(inp, 0, -1)

        poles = torch.as_tensor(pole, dtype=inp.dtype, device=inp.device)
        poles = poles.pow(torch.arange(1, n-1, dtype=inp.dtype, device=inp.device))
        poles = poles + (polen * polen) / poles

        out = torch.matmul(inp.unsqueeze(-2), poles.unsqueeze(-1)).squeeze(-1)
        out = out + inp0.unsqueeze(-1)
        if keepdim:
            out = movedim1(out, -1, dim)
        else:
            out = out.squeeze(-1)

        pole = pole ** (max_iter - 1)
        out = out / (1 - pole * pole)

    return out


@torch.jit.script
def dct2_initial(inp, pole: float, dim: int = -1, keepdim: bool = False):
    # Ported from scipy:
    # https://github.com/scipy/scipy/blob/master/scipy/ndimage/src/ni_splines.c
    #
    # I (YB) unwarped and simplied the terms so that I could use a dot
    # product instead of a loop.
    # It should certainly be possible to derive a version for max_iter < n,
    # as JA did for DCT1, to avoid long recursions when `n` is large. But
    # I think it would require a more complicated anticausal/final condition.

    n = inp.shape[dim]

    polen = pole ** n
    pole_last = polen * (1 + 1/(pole + polen * polen))
    inp00 = inp[0]
    inp0 = inp[0] + pole_last * inp[-1]
    inp = inp[1:-1]
    inp = movedim1(inp, 0, -1)

    poles = torch.as_tensor(pole, dtype=inp.dtype, device=inp.device)
    poles = (poles.pow(torch.arange(1, n-1, dtype=inp.dtype, device=inp.device)) +
             poles.pow(torch.arange(2*n-2, n, -1, dtype=inp.dtype, device=inp.device)))

    out = torch.matmul(inp.unsqueeze(-2), poles.unsqueeze(-1)).squeeze(-1)

    out = out + inp0.unsqueeze(-1)
    out = out * (pole / (1 - polen * polen))
    out = out + inp00.unsqueeze(-1)

    if keepdim:
        out = movedim1(out, -1, dim)
    else:
        out = out.squeeze(-1)

    return out


@torch.jit.script
def dft_final(inp, pole: float, dim: int = -1, keepdim: bool = False):

    assert inp.shape[dim] > 1
    max_iter: int = int(math.ceil(-30./math.log(abs(pole))))
    max_iter = min(max_iter, inp.shape[dim])

    poles = torch.as_tensor(pole, dtype=inp.dtype, device=inp.device)
    poles = poles.pow(torch.arange(2, max_iter+1, dtype=inp.dtype, device=inp.device))

    inp = movedim1(inp, dim, 0)
    inp0 = inp[-1]
    inp = inp[:max_iter-1]
    inp = movedim1(inp, 0, -1)
    out = torch.matmul(inp.unsqueeze(-2), poles.unsqueeze(-1)).squeeze(-1)
    out = out.add(inp0.unsqueeze(-1), alpha=pole)
    if keepdim:
        out = movedim1(out, -1, dim)
    else:
        out = out.squeeze(-1)

    pole = pole ** max_iter
    out = out / (pole - 1)
    return out


@torch.jit.script
def dct1_final(inp, pole: float, dim: int = -1, keepdim: bool = False):
    inp = movedim1(inp, dim, 0)
    out = pole * inp[-2] + inp[-1]
    out = out * (pole / (pole*pole - 1))
    if keepdim:
        out = movedim1(out.unsqueeze(0), 0, dim)
    return out


@torch.jit.script
def dct2_final(inp, pole: float, dim: int = -1, keepdim: bool = False):
    # Ported from scipy:
    # https://github.com/scipy/scipy/blob/master/scipy/ndimage/src/ni_splines.c
    inp = movedim1(inp, dim, 0)
    out = inp[-1] * (pole / (pole - 1))
    if keepdim:
        out = movedim1(out.unsqueeze(0), 0, dim)
    return out


@torch.jit.script
class CoeffBound:

    def __init__(self, bound: int):
        self.bound = bound

    def initial(self, inp, pole: float, dim: int = -1, keepdim: bool = False):
        if self.bound in (0, 2):    # zero, dct1
            return dct1_initial(inp, pole, dim, keepdim)
        elif self.bound in (1, 3):  # nearest, dct2
            return dct2_initial(inp, pole, dim, keepdim)
        elif self.bound == 6:       # dft
            return dft_initial(inp, pole, dim, keepdim)
        else:
            raise NotImplementedError

    def final(self, inp, pole: float, dim: int = -1, keepdim: bool = False):
        if self.bound in (0, 2):    # zero, dct1
            return dct1_final(inp, pole, dim, keepdim)
        elif self.bound in (1, 3):  # nearest, dct2
            return dct2_final(inp, pole, dim, keepdim)
        elif self.bound == 6:       # dft
            return dft_final(inp, pole, dim, keepdim)
        else:
            raise NotImplementedError


@torch.jit.script
def filter(inp, bound: CoeffBound, poles: List[float],
           dim: int = -1, inplace: bool = False):

    if not inplace:
        inp = inp.clone()

    if inp.shape[dim] == 1:
        return inp

    gain = get_gain(poles)
    inp *= gain
    inp = movedim1(inp, dim, 0)
    n = inp.shape[0]

    for pole in poles:
        inp[0] = bound.initial(inp, pole, dim=0, keepdim=False)

        for i in range(1, n):
            inp[i].add_(inp[i-1], alpha=pole)

        inp[-1] = bound.final(inp, pole, dim=0, keepdim=False)

        for i in range(n-2, -1, -1):
            inp[i].neg_().add_(inp[i+1]).mul_(pole)

    inp = movedim1(inp, 0, dim)
    return inp


@torch.jit.script
def spline_coeff(inp, bound: int, order: int, dim: int = -1,
                 inplace: bool = False):
    """Compute the interpolating spline coefficients, for a given spline order
    and boundary conditions, along a single dimension.

    Parameters
    ----------
    inp : tensor
    bound : {2: dct1, 6: dft}
    order : {0..7}
    dim : int, default=-1
    inplace : bool, default=False

    Returns
    -------
    out : tensor

    """
    if not inplace:
        inp = inp.clone()

    if order in (0, 1):
        return inp

    poles = get_poles(order)
    return filter(inp, CoeffBound(bound), poles, dim=dim, inplace=True)


@torch.jit.script
def spline_coeff_nd(inp, bound: List[int], order: List[int],
                    dim: Optional[int] = None, inplace: bool = False):
    """Compute the interpolating spline coefficients, for a given spline order
    and boundary condition, along the last `dim` dimensions.

    Parameters
    ----------
    inp : (..., *spatial) tensor
    bound : List[{2: dct1, 6: dft}]
    order : List[{0..7}]
    dim : int, default=`inp.dim()`
    inplace : bool, default=False

    Returns
    -------
    out : (..., *spatial) tensor

    """
    if not inplace:
        inp = inp.clone()

    if dim is None:
        dim = inp.dim()

    bound = pad_list_int(bound, dim)
    order = pad_list_int(order, dim)

    for d, b, o in zip(range(dim), bound, order):
        inp = spline_coeff(inp, b, o, dim=-dim + d, inplace=True)

    return inp
