# first import basic IO tools
from . import term
from .logging import *
from .parser import ArgParser

# primary utilities
from .utility import *

# c bindings
try:
    from . import bindings
except ImportError:
    # check if importing from a local repository - most likely the bindings library
    # has not been built or an incorrect python version is used, so print out
    # a message to hopefully get the developer on the right track
    import os.path as path
    if path.exists(path.abspath(path.join(__file__ ,"../../CMakeLists.txt"))):
        error('cannot import freesurfer C bindings. Make sure that the `python/bindings` '
              'subdirectory has been built and the current python version matches the '
              'one used to compile the bindings library')
    raise

from .transform import *

from .surface import Surface
from .volume import Volume
from .bindings import read_annotation

from .geometry import *
from ._surface import *
from ._normalize import *
from .freeview import *
from .mrisp import *
from . import metrics

