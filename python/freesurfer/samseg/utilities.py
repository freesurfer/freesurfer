import json
import datetime
import numpy as np
from _operator import itemgetter
from collections import namedtuple


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


def kvlReadSharedGMMParameters(fileName):
    # Read text file where each line is formatted as:
    #   mergedName numberOfComponents searchString(s)
    print('Reading contexts of %s' % fileName)
    sharedGMMParameters = []
    with open(fileName) as fid:
        for textLine in fid.readlines():
            # Remove leading and trailing blanks
            components = textLine.strip().split()
            if len(components) > 0 and components[0] != '#':
                if len(components) < 2:
                    raise ValueError( 'Can''t parse line: {0}'.format(textLine))
                mergedName = components[0]
                numberOfComponents = int(components[1])
                searchStrings = components[2:]
                # Add to sharedGMMParameters structure array
                GMMparameter = namedtuple('GMMparameter', 'mergedName numberOfComponents searchStrings')
                sharedGMMParameters.append(GMMparameter(mergedName, numberOfComponents, searchStrings))
    return sharedGMMParameters


def kvlReadCompressionLookupTable(fileName):
    # Format is "FreeSurferLabel compressedLabel name R G B A"
    print('Reading contexts of %s' % fileName)
    table = []
    with open(fileName) as fid:
        for line in fid.readlines():
            FreeSurferLabel, compressedLabel, name, R, G, B, A = [
                data_type(value) for data_type, value in zip(
                    (int, int, str, int, int, int, int),
                    line.split())]
            # Add contents to output matrices
            table.append({
                'FreeSurferLabel': FreeSurferLabel,
                'compressedLabel': compressedLabel,
                'name': name,
                'color': [R, G, B, A],
            })
    # Sort output according to compressedLabel
    table = sorted(table, key=itemgetter('compressedLabel'))
    FreeSurferLabels, names, colors = [[entry[key] for entry in table] for key in ['FreeSurferLabel', 'name', 'color']]
    return FreeSurferLabels, names, colors
