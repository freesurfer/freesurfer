import numpy as np
from . import bindings, LinearTransform


class Surface(bindings.Surface):

    def surf2vox(self, vol):
        return LinearTransform(self._compute_surf2vox(vol))

    def vox2surf(self, vol):
        return self.surf2vox(vol).inverse()
