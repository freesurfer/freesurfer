try:
    from .fsbindings import *
except ImportError:
    raise ImportError('fsbindings have not been built for this python version')
