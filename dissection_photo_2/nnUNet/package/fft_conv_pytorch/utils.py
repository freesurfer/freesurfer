from timeit import Timer
from typing import Callable, NamedTuple

import numpy as np
import torch
from torch import Tensor


class Benchmark(NamedTuple):
    mean: float
    std: float

    def __repr__(self):
        return f"BenchmarkResult(mean: {self.mean:.3e}, std: {self.std:.3e})"

    def __str__(self):
        return f"({self.mean:.3e} \u00B1 {self.std:.3e}) s"


def benchmark(fn: Callable, *args, num_iterations: int = 10, **kwargs) -> Benchmark:
    timer = Timer(
        "fn(*args, **kwargs)",
        globals={"fn": fn, "args": args, "kwargs": kwargs},
    )
    times = timer.repeat(number=1, repeat=num_iterations + 1)
    return Benchmark(np.mean(times[1:]).item(), np.std(times[1:]).item())


def _assert_almost_equal(x: Tensor, y: Tensor) -> bool:
    abs_error = torch.abs(x - y)
    assert abs_error.mean().item() < 5e-5
    assert abs_error.max().item() < 1e-4
    return True


def _gcd(x: int, y: int) -> int:
    while y:
        x, y = y, x % y
    return x
