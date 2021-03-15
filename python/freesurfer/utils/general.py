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
        # NOTE: The reason we specify h5py below is to deal with a bug in
        # this version of TF. It shouldn't be needed with newer TF versions
        error('This tool requires TensorFlow, but fspython does not ship with tensorflow')
        print('by default. You (or a sys admin) can install the CPU version via:\n')
        print('  fspython -m pip install h5py==2.10.0 tensorflow==1.13.1\n')
        print('or the GPU version via:\n')
        print('  fspython -m pip install h5py==2.10.0 tensorflow-gpu==1.13.1\n')
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
