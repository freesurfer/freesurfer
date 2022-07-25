import numpy as np

from . import lookups
from .lookups import LookupTable
from .lookups import RecodingLookupTable
from .ndarray import ArrayContainerTemplate
from .deprecations import deprecate, replace, unsure, notneeded


@replace('sf.labels.recode')
def recode(seg, mapping):
    """
    Recodes a hard segmentation via a label mapping.

    Parameters:
        seg: Segmentation to recode.
        mapping: Label index map. Can be list, dict, or RecodingLookupTable.
    Returns:
        Recoded array.
    """

    # this is such an ugly hack - we really shouldn't include
    # this kind of support 
    #if seg.__class__.__name__ in ('Tensor', 'EagerTensor'):
    #    import neurite as ne
    #    return ne.utils.seg.recode(seg, mapping)

    seg_data = seg.data if isinstance(seg, ArrayContainerTemplate) else seg
    recoded = np.zeros_like(seg_data, dtype=np.int32)
    recoded_lut = None

    # convert mapping to dictionary
    if isinstance(mapping, np.ndarray) and mapping.ndim != 1:
        raise ValueError('Label recoding list must be 1D, but got %dD array.' % mapping.ndim)
    elif isinstance(mapping, (list, tuple, np.ndarray)):
        mapping = {l: i + 1 for i, l in enumerate(mapping)}
    elif isinstance(mapping, RecodingLookupTable):
        recoded_lut = mapping.target_lut
        mapping = mapping.mapping
    elif not isinstance(mapping, dict):
        raise ValueError('Invalid mapping type \'%s\'.' % type(mapping).__name__)

    # do the recoding
    for i in np.unique(seg_data):
        new_i = mapping.get(i)
        if new_i is not None:
            # TODO: this is probably not the most efficient method
            recoded[seg_data == i] = new_i

    # convert back
    if isinstance(seg, ArrayContainerTemplate):
        recoded = seg.copy(recoded)
        recoded.lut = recoded_lut

    return recoded


@replace('seg.onehot()')
def hard_to_prob(seg, labels, dtype='float32'):
    """
    Converts a hard segmentation to a probabilistic segmentation.
    """
    seg_data = seg.data if isinstance(seg, ArrayContainerTemplate) else seg

    if isinstance(labels, int):
        labels = np.arange(labels) + 1

    # make prob seg
    prob = np.zeros((*seg_data.shape, len(labels) + 1), dtype=dtype)
    for i, l in enumerate(labels):
        prob[..., i + 1] = seg_data == l
    prob[..., 0] = prob[..., 1:].sum(axis=-1) == 0

    # convert back
    if isinstance(seg, ArrayContainerTemplate):
        prob = seg.copy(prob)
        prob.lut = None

    return prob


@replace('seg.collapse()')
def prob_to_hard(seg, labels=None):
    """
    Converts a probabilistic segmentation to a hard segmentation.
    """
    seg_data = seg.data if isinstance(seg, ArrayContainerTemplate) else seg

    hard = np.argmax(seg_data, axis=-1)
    if labels is not None:
        if len(labels) == seg_data.shape[-1] - 1:
            mapping = {i + 1: l for i, l in enumerate(labels)}
        elif len(labels) == seg_data.shape[-1]:
            mapping = {i: l for i, l in enumerate(labels)}
        else:
            raise ValueError('Label list does not match prob seg frames.')
        hard = recode(hard, mapping)

    # convert back
    if isinstance(seg, ArrayContainerTemplate):
        hard = seg.copy(hard)

    return hard


@replace('sf.labels.dice')
def dice(a, b, labels=None):
    '''
    Computes dice coefficients for each label between two hard segmentations,
    using the formula:

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
    Computes jaccard coefficients for each label between two hard segmentations,
    using the formula:
    
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
