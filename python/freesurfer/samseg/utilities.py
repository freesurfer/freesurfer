import json
import numpy as np


def forceFortranOrder(funcname):
    old_func = getattr(np, funcname)
    def new_func(*args, **kwargs):
        kwargs['order'] = 'F'
        return old_func(*args, **kwargs)
    setattr(np, funcname, new_func)

for funcname in ('array', 'zeros', 'empty', 'zeros_like', 'empty_like'):
    forceFortranOrder(funcname)


def dump_dict(value):
    return getattr(value, 'dump_dict', getattr(value, '__dict__', value))


class Specification:
    def __init__(self, params):
        for key, value in params.items():
            self.__setattr__(key, value)

    @property
    def dump_dict(self):
        return {key: dump_dict(value) for key, value in self.__dict__.items()}

    def toJSON(self):
        return json.dumps(self, default=dump_dict, sort_keys=True)


def requireNumpyArray(np_array):
    return np.require(np_array, requirements=['F_CONTIGUOUS', 'ALIGNED'])


def ensureDims(np_array, dims):
    if np_array.ndim < dims:
        return ensureDims(np.expand_dims(np_array, axis=dims), dims)
    elif np_array.ndim == dims:
        return np_array


def icv(structures, includeStructures=None):
    if not includeStructures:
        print("using default intracranial structures to compute sbtiv measure")
        includeStructures = [  
            'Brain-Stem',
            'CSF',
            'Left-Cerebellum-Cortex',
            'Right-Cerebellum-Cortex',
            '4th-Ventricle',
            'Left-Cerebral-Cortex',
            'Right-Cerebral-Cortex',
            'Left-Cerebellum-White-Matter',
            'Right-Cerebellum-White-Matter',
            'Right-Cerebral-White-Matter',
            'Left-Cerebral-White-Matter',
            'Right-choroid-plexus',
            'Right-Amygdala',
            'Left-Amygdala',
            'Right-Hippocampus',
            'Right-Inf-Lat-Vent',
            'Left-Hippocampus',
            'Left-choroid-plexus',
            'Left-Inf-Lat-Vent',
            'Optic-Chiasm',
            'Right-VentralDC',
            'Left-VentralDC',
            '3rd-Ventricle',
            'Left-Accumbens-area',
            'Left-Putamen',
            'Right-Putamen',
            'Right-vessel',
            'Right-Accumbens-area',
            'WM-hypointensities',
            'Left-Lateral-Ventricle',
            'Left-vessel',
            'Right-Lateral-Ventricle',
            'Right-Pallidum',
            'Left-Caudate',
            'Right-Thalamus-Proper',
            'Left-Pallidum',
            'Right-Caudate',
            'Left-Thalamus-Proper',
            'non-WM-hypointensities',
            '5th-Ventricle'
        ]
    return sum(structure[1] for structure in structures if structure[0] in includeStructures)