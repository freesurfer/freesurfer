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
