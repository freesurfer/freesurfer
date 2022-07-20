import os
import sys

# this is where we keep the archived version
archived = '/autofs/space/freesurfer/python/fsmodule'
if not os.path.exists(archived):
    raise ImportError('the freesurfer package has been deprecated in favor of surfa')

sys.path.insert(0, archived)

from fsarchive import *

import warnings

def formatwarning(message, category, filename, lineno, line=None):
    return f'{category.__name__}: {message}\n'
warnings.formatwarning = formatwarning

warnings.warn('the freesurfer package has been deprecated in favor of surfa, but loading the archived version')
