import importlib.util
spec = importlib.util.find_spec('fsbindings.fsbindings')
if (spec):
    from .fsbindings import *
#else:
#    import warnings
#    warnings.warn('fsbindings c++ libraries have not been built for this python version')
