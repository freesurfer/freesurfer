import os
import scipy.ndimage
import numpy as np
import surfa as sf

from freesurfer import samseg


def run(cmd):
    """
    Run a command in a bash shell. Output is silenced, but is printed if an error occurs.
    """
    print(f'Running command: {cmd}')
    output, ret = sf.system.collect_output(cmd)
    if ret != 0:
        print(output)
        sf.system.fatal('Command failed')


def spherical_strel(radius, pixsize=1.0):
    """
    Compute a 3D spherical binary structure for mask manipulation.
    """
    pixsize = np.array([pixsize] * 3)
    shape = np.ceil(2 * radius / pixsize + 1).astype(int)
    shape += np.mod(shape + 1, 2)
    center = (shape - 1) / 2
    coords = np.array(np.ogrid[:shape[0], :shape[1], :shape[2]], dtype=object)
    return np.sum((coords - center) ** 2, axis=0) <= (radius ** 2)


def read_compression_lookup_table(filename):
    """
    Read a compressed label lookup table file into an ordered dictionary
    to labels and names. This also returns corresponding label indices and names
    in a tuple, although we can probably re-extract this info from the labelMapping
    object down the road.
    """
    labelMapping = sf.LabelLookup()
    labels, names, colors = samseg.kvlReadCompressionLookupTable(filename)
    labels = np.array(labels)
    for label, name, color in zip(labels, names, colors):
        labelMapping[label] = (name, color)
    return (labelMapping, names, labels)


def get_largest_cc(mask):
    """
    Find the largest connected component of a binary mask. All over components are
    masked away in the returned array.
    ATH TODO: This should be implemented as a function of the Volume object.
    """
    labels = scipy.ndimage.label(mask)[0]
    return labels == np.argmax(np.bincount(labels.flatten())[1:]) + 1


def geometries_differ(a, b):
    """
    Compare the similarity of two volume geometries.
    """
    if not np.array_equal(a.shape[:3], b.shape[:3]):
        return True
    if not np.max(np.abs(a.voxsize - b.voxsize)) > 1e-5:
        return True
    if not np.max(np.abs(a.matrix - b.matrix)) > 1e-5:
        return True
    return False
