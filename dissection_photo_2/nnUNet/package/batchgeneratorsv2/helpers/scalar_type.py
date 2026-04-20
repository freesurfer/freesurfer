from typing import Union, Tuple, Callable
import numpy as np


RandomScalar = Union[int, float, Tuple[float, float], Callable[..., Union[int, float]]]


def sample_scalar(scalar_type: RandomScalar, *args, **kwargs):
    if isinstance(scalar_type, (int, float)):
        return scalar_type
    elif isinstance(scalar_type, (list, tuple)):
        assert len(scalar_type) == 2, 'if list is provided, its length must be 2'
        assert scalar_type[0] <= scalar_type[1], 'if list is provided, first entry must be smaller or equal than second entry, ' \
                                                'otherwise we cannot sample using np.random.uniform'
        if scalar_type[0] == scalar_type[1]:
            return scalar_type[0]
        return np.random.uniform(*scalar_type)
    elif callable(scalar_type):
        return scalar_type(*args, **kwargs)
    else:
        raise RuntimeError('Unknown type: %s. Expected: int, float, list, tuple, callable', type(scalar_type))


if __name__ == '__main__':
    sample_scalar(0.5)
    sample_scalar((0, 1))
    sample_scalar(lambda: np.random.uniform(-1, 2))
    sample_scalar(lambda x, y: np.random.uniform(x, y), 0.5, 2)