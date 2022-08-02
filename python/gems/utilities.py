import json
import numpy as np
from copy import deepcopy


def forceFortranOrder(funcname):
    old_func = getattr(np, funcname)
    def new_func(*args, **kwargs):
        kwargs['order'] = 'F'
        return old_func(*args, **kwargs)
    setattr(np, funcname, new_func)


class Specification:
    def __init__(self, params):
        self.depth = 1
        for key, value in params.items():
            self.__setattr__(key, value)

    def __str__(self):
        spec = '{\n'
        for key, value in self.__dict__.items():
            spec += '%s%s: ' % (self.depth * '    ', key)
            if isinstance(value, Specification):
                value.depth = self.depth + 1
                spec += self._indented(value)
            if isinstance(value, list):
                spec += ', '.join([self._indented(i) for i in value]) + '\n'
            else:
                spec += str(value) + '\n'
            if key == '':
                print(type(variable_name))
        return '%s%s}' % (spec, (self.depth - 1) * '    ')

    def _indented(self, item):
        if isinstance(item, Specification):
            item.depth = self.depth + 1
            return str(item)
        else:
            return str(item)

    def merged(self, params):
        newSpec = deepcopy(self)
        if not params:
            return newSpec
        for key, value in params.items():
            if not hasattr(newSpec, key):
                raise ValueError("specification does not contain parameter '%s'" % key)
            newSpec.__setattr__(key, value)
        return newSpec


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
            'Right-Thalamus',
            'Left-Pallidum',
            'Right-Caudate',
            'Left-Thalamus',
            'non-WM-hypointensities',
            '5th-Ventricle',
            'Lesions'
        ]
    return sum(structure[1] for structure in structures if structure[0] in includeStructures)
