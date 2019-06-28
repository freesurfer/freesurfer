import numpy as np
from . import bindings, LinearTransform


class Volume(bindings.Volume):

    def vox2ras(self):
        return LinearTransform(self.affine)

    def ras2vox(self):
        return self.vox2ras().inverse()

    def vox2surf(self):
        return LinearTransform(self._compute_vox2surf())

    def surf2vox(self):
        return self.vox2surf().inverse()
