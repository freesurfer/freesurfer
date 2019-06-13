import numpy as np
from . import bindings, LinearTransform


class Surface(bindings.Surface):

    def surf2vox(self, vol):
        return LinearTransform(self._compute_surf2vox(vol))

    def vox2surf(self, vol):
        return self.surf2vox(vol).inverse()


def sample_parameterization(sp, mapping):
    mapping = np.moveaxis(mapping, -1, 0)
    weights = mapping[0]
    ui = mapping[1].astype(int)
    vi = mapping[2].astype(int)
    return np.sum(weights * sp[ui, vi], axis=1)
