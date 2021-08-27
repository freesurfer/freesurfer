# -- freesurfer --
# neuroimaging analysis and visualization suite

# import bindings used to wrap c++ routines
# if this import fails, it probably means that either the library has
# not been built or you're python version is different from what was
# originally used to compile the bindings
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

# TODO: these should be deprecated
from .normalize import *
from .deprecations import *
