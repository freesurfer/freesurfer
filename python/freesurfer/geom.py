import numpy as np
import scipy.ndimage
from collections.abc import Iterable


class Slicing(tuple):
    '''Slice tuple used for indexing subregions of a numpy array.'''

    def __new__(cls, start, stop):
        if len(start) != len(stop):
            raise ValueError('start coord (%dD) does not match stop coord (%dD))' % (len(start), len(stop)))
        return super(Slicing, cls).__new__(cls, [slice(s, t) for s, t in zip(start, stop)])

    @property
    def start(self):
        return tuple([s.start for s in self])

    @property
    def stop(self):
        return tuple([s.stop for s in self])

    @property
    def shape(self):
        return tuple([s.stop - s.start for s in self])

    def grow(self, dist):
        if not isinstance(dist, Iterable):
            dist = [dist] * len(self)
        start = [s.start - int(d) for s, d in zip(self, dist)]
        stop  = [s.stop  + int(d) for s, d in zip(self, dist)]
        return Slicing(start, stop)

    def shrink(self, dist):
        if isinstance(dist, Iterable):
            dist = np.array(dist)
        return self.grow(dist * -1)


def bbox(mask):
    '''Bounding box around the object in a binary image.'''
    return scipy.ndimage.find_objects(mask)[0]


def cmass(image):
    '''Center of mass of an image.'''
    return scipy.ndimage.center_of_mass(image)
