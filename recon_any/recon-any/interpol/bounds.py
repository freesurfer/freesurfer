import torch
from enum import Enum
from typing import Optional
from .jit_utils import floor_div
Tensor = torch.Tensor


class BoundType(Enum):
    zero = zeros = 0
    replicate = nearest = 1
    dct1 = mirror = 2
    dct2 = reflect = 3
    dst1 = antimirror = 4
    dst2 = antireflect = 5
    dft = wrap = 6


class ExtrapolateType(Enum):
    no = 0     # threshold: (0, n-1)
    yes = 1
    hist = 2   # threshold: (-0.5, n-0.5)


@torch.jit.script
class Bound:

    def __init__(self, bound_type: int = 3):
        self.type = bound_type

    def index(self, i, n: int):
        if self.type in (0, 1):  # zero / replicate
            return i.clamp(min=0, max=n-1)
        elif self.type in (3, 5):  # dct2 / dst2
            n2 = n * 2
            i = torch.where(i < 0, (-i-1).remainder(n2).neg().add(n2 - 1),
                            i.remainder(n2))
            i = torch.where(i >= n, -i + (n2 - 1), i)
            return i
        elif self.type == 2:  # dct1
            if n == 1:
                return torch.zeros(i.shape, dtype=i.dtype, device=i.device)
            else:
                n2 = (n - 1) * 2
                i = i.abs().remainder(n2)
                i = torch.where(i >= n, -i + n2, i)
                return i
        elif self.type == 4:  # dst1
            n2 = 2 * (n + 1)
            first = torch.zeros([1], dtype=i.dtype, device=i.device)
            last = torch.full([1], n - 1, dtype=i.dtype, device=i.device)
            i = torch.where(i < 0, -i - 2, i)
            i = i.remainder(n2)
            i = torch.where(i > n, -i + (n2 - 2), i)
            i = torch.where(i == -1, first, i)
            i = torch.where(i == n, last, i)
            return i
        elif self.type == 6:  # dft
            return i.remainder(n)
        else:
            return i

    def transform(self, i, n: int) -> Optional[Tensor]:
        if self.type == 4:  # dst1
            if n == 1:
                return None
            one = torch.ones([1], dtype=torch.int8, device=i.device)
            zero = torch.zeros([1], dtype=torch.int8, device=i.device)
            n2 = 2 * (n + 1)
            i = torch.where(i < 0, -i + (n-1), i)
            i = i.remainder(n2)
            x = torch.where(i == 0, zero, one)
            x = torch.where(i.remainder(n + 1) == n, zero, x)
            i = floor_div(i, n+1)
            x = torch.where(torch.remainder(i, 2) > 0, -x, x)
            return x
        elif self.type == 5:  # dst2
            i = torch.where(i < 0, n - 1 - i, i)
            x = torch.ones([1], dtype=torch.int8, device=i.device)
            i = floor_div(i, n)
            x = torch.where(torch.remainder(i, 2) > 0, -x, x)
            return x
        elif self.type == 0:  # zero
            one = torch.ones([1], dtype=torch.int8, device=i.device)
            zero = torch.zeros([1], dtype=torch.int8, device=i.device)
            outbounds = ((i < 0) | (i >= n))
            x = torch.where(outbounds, zero, one)
            return x
        else:
            return None
