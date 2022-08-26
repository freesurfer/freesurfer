# This is the original freesurfer python package, which has no been
# deprecated in favor of the surfa library. While we no longer distribute
# this library in freesurfer distributions, we would like to continue to
# make it available so that things don't start to break. This folder should
# be completely removed at some point. If you're reading this past 2022, talk
# to Bruce about deleting it. This package depends on c++ bindings that have since
# been converted to brigde surfa objects (see `fsbindings` folder). It is no
# longer possible to build these original bindings libraries, unless you checkout
# the commit f73b799c and run `make` in the `python/bindings` subdirectory. So,
# to allow for the original bindings to be imported, we've added a few prebuilt
# `bindings.cpython*` libraries to the repository. Note that these are only valid
# for linux, and will not be importable on macos.

# Let's only allow import on the martinos network, otherwise throw an error.
import os
if not os.path.exists('/autofs/space/freesurfer/python'):
    raise ImportError('The freesurfer package has been deprecated in favor of surfa')

import warnings

def formatwarning(message, category, filename, lineno, line=None):
    return f'{category.__name__}: {message}\n'
warnings.formatwarning = formatwarning

# setting the FS_SURFA_PORT will print out suggested function for porting to the surfa library
warnings.warn('freesurfer package: set FS_SURFA_PORT env var to print surfa porting suggestions')

from . import bindings

# general utilities
from . import utils
from .utils import run
from .utils import collect_output
from .utils import warning
from .utils import error
from .utils import fatal
from .utils import fshome

# geometric utilities and algorithms
from . import metrics
from . import geom
from .transform import Geometry
from .transform import LinearTransform
from .transform import Warp

# label
from . import label
from .label import LookupTable

# ND array containers
from .ndarray import Overlay
from .ndarray import Image
from .ndarray import Volume

# surfaces
from .surface import Surface

# visualization
from .freeview import fv
from .freeview import fvoverlay
from .freeview import Freeview

from .normalize import *
from .deprecations import *
