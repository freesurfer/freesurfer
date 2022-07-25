import numpy as np
from scipy.ndimage.morphology import binary_erosion
from scipy.spatial.distance import cdist

from .deprecations import deprecate, replace, unsure, notneeded, notimplemented


@replace('sf.labels.dice')
def dice(a, b, labels=None):
    '''
    Computes dice coefficients for each label between two hard segmentations.
    
        dice = 2 * |A ∩ B| / |A| + |B|

    If labels are not provided, all unique integers are used except for 0.
    '''
    if labels is None:
        labels = np.unique(np.concatenate([a, b]))
        labels = np.delete(labels, np.where(labels == 0))
    
    result = {}
    for l in labels:
        mask1 = a == l
        mask2 = b == l
        top = 2.0 * np.logical_and(mask1, mask2).sum()
        bottom = mask1.sum() + mask2.sum()
        bottom = np.maximum(bottom, np.finfo(float).eps)  # add epsilon
        result[l] = top / bottom
    
    return result


@replace('sf.labels.jaccard')
def jaccard(a, b, labels=None):
    '''
    Computes jaccard coefficients for each label between two hard segmentations.
    
        jaccard = |A ∩ B| / |A ∪ B|

    If labels are not provided, all unique integers are used except for 0.
    '''
    if labels is None:
        labels = np.unique(np.concatenate([a, b]))
        labels = np.delete(labels, np.where(labels == 0))
    
    result = {}
    for l in labels:
        mask1 = a == l
        mask2 = b == l
        top = np.logical_and(mask1, mask2).sum()
        bottom = np.logical_or(mask1, mask2).sum()
        bottom = np.maximum(bottom, np.finfo(float).eps)  # add epsilon
        result[l] = top / bottom
    
    return result


@notimplemented
def hausdorff(a, b, thresh=0.01, measure=np.max):
    '''
    Hausdorff distance between two volumes. Computes the max of
    the minimum euclidian distances between bordering voxels of
    the two input volumes.

    Args:
        thresh: Input threshold for masking. Defaults to 0.01.
        measure: Measure function of reported distance. Defaults to `np.max`.
    '''
    mask1 = a >= thresh
    mask2 = b >= thresh
    
    # find borders by eroding binary masks
    border1 = np.logical_and(mask1, np.logical_not(binary_erosion(mask1)))
    border2 = np.logical_and(mask2, np.logical_not(binary_erosion(mask2)))
    
    # compute voxel distances
    dist = cdist(np.vstack(np.nonzero(border1)).T, np.vstack(np.nonzero(border2)).T)
    mindist = np.concatenate((np.amin(dist, axis=0), np.amin(dist, axis=1)))
    return measure(mindist)
