from .log import *
from .util import *

# these require 3rd party dependencies
# make sure they're installed before importing
try:
    from . import metrics
except ModuleNotFoundError as e:
    warning(e.msg + ", used in freesurfer python package 'metrics'")
try:
    from .freeview import fv
except ModuleNotFoundError as e:
    warning(e.msg + ", used in freesurfer python package 'freeview'")
try:
    from ._normalize import *
except ModuleNotFoundError as e:
    warning(e.msg + ", used in freesurfer python package '_normalize'")

try:
    from ._surface import *
except ModuleNotFoundError as e:
    warning(e.msg + ", used in freesurfer python package '_surface'")
