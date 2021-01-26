import os
import sys
import pickle
import numpy as np
import datetime as dt
import traceback

from .logging import error


def fshome():
    '''Returns the freesurfer home directory.'''
    return os.environ.get('FREESURFER_HOME')


def readlines(filename):
    '''Reads the lines of a text file in to a list.'''
    with open(filename) as file:
        content = file.read().splitlines()
    return content


def write_pickle(item, filename):
    '''Pickles item (with highest protocol) to a file.'''
    with open(filename, 'wb') as file:
        pickle.dump(item, file, protocol=pickle.HIGHEST_PROTOCOL)


def read_pickle(filename):
    '''Reads pickled item from file.'''
    with open(filename, 'rb') as file:
        item = pickle.load(file)
    return item


def check_tensorflow():
    '''Ensures that tensorflow is installed in fspython.'''
    try:
        import tensorflow
    except ModuleNotFoundError:
        error('This tool requires TensorFlow, but fspython does not ship with tensorflow')
        print('by default. You (or a sys admin) can install the CPU version via:\n')
        print('  fspython -m pip install tensorflow==1.13.1\n')
        print('or the GPU version via:\n')
        print('  fspython -m pip install tensorflow-gpu==1.13.1\n')
        sys.exit(1)
    except ImportError as err:
        print(traceback.format_exc())
        error('Found fspython TensorFlow, but could not import (see error above).')
        print('If you\'ve install the GPU version, please ensure that the appropriate')
        print('CUDA libraries are installed on your system and available in your environment.')
        sys.exit(1)


class Timer:
    '''A simple timer class to track process speed.'''
    def __init__(self, message=None):
        if message: print(message)
        self.start_time = dt.datetime.now()

    @property
    def elapsed(self):
        return dt.datetime.now() - self.start_time

    def mark(self, message):
        print('%s: %s' % (message, str(self.elapsed)))


class LookupTable(dict):
    # TODOC

    class Element:
        def __init__(self, name, color):
            self.name = name
            if color is None:
                color = np.random.randint(0, 256, size=(3))
            if len(color) == 3:
                color = np.append(color, 255)
            elif len(color) != 4:
                raise ValueError('Color must be a 4-element RGBA uchar array')
            self.color = np.array(color, dtype='uint8')

    def __str__(self):
        col1 = len(str(max(self.keys()))) + 1
        col2 = max([len(elt.name) for elt in self.values()]) + 2
        lines = []
        for idx, elt in self.items():
            colorstr = '(' + ' '.join([str(c).ljust(3) for c in elt.color]) + ')'
            lines.append(str(idx).ljust(col1) + elt.name.ljust(col2) + colorstr)
        return '\n'.join(lines)

    def add(self, index, name, color):
        self[index] = LookupTable.Element(name, color)

    def search(self, name):
        allcaps = name.upper()
        return [idx for idx, elt in self.items() if allcaps in elt.name.upper()]

    @classmethod
    def read(cls, filename):
        lut = cls()
        with open(filename, 'r') as file:
            lines = file.readlines()
        for line in lines:
            split = line.lstrip().split()
            if split and not split[0].startswith('#'):
                idx, name = split[:2]
                if len(split) >= 5:
                    color = list(map(int, split[2:6]))
                    color[3] = 255 - color[3]  # invert alpha value
                else:
                    color = None
                lut.add(int(idx), name, color)
        return lut

    @classmethod
    def read_default(cls):
        return cls.read(os.path.join(fshome(), 'FreeSurferColorLUT.txt'))

    def write(self, filename):
        col1 = len(str(max(self.keys()))) + 1  # find largest index
        col2 = max([len(elt.name) for elt in self.values()]) + 2  # find longest name
        with open(filename, 'w') as file:
            file.write('#'.ljust(col1) + 'Label Name'.ljust(col2) + 'R   G   B   A\n\n')
            for idx, elt in self.items():
                color = elt.color
                color[3] = 255 - color[3]  # invert alpha value
                colorstr = ' '.join([str(c).ljust(3) for c in color])
                file.write(str(idx).ljust(col1) + elt.name.ljust(col2) + colorstr + '\n')
