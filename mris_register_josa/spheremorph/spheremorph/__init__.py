# ---- spheremorph ----
# unsupervised learning for image registration on spheres

__version__ = '0.1'

import voxelmorph

# import backend-dependent submodules
backend = voxelmorph.py.utils.get_backend()

if backend == 'pytorch':
    pass
    # from . import torch
    # from .torch import layers
    # from .torch import networks
    # from .torch import losses
else:
    from . import tf
    from .tf import layers
    from .tf import networks
    from .tf import losses
    from .tf import utils
    from . import py
