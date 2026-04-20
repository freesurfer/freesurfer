from typing import Union, Tuple, List, Callable

import numpy as np
from skimage.measure import label
from skimage.morphology import ball
from skimage.transform import resize
#import cc3d


def generate_ball(radius: Union[Tuple, List], spacing: Union[Tuple, List] = (1, 1, 1), dtype=np.uint8) -> np.ndarray:
    """
    Returns a ball/ellipsoid corresponding to the specified size (radius = list/tuple of len 3 with one radius per axis)
    If you use spacing, both radius and spacing will be interpreted relative to each other, so a radius of 10 with a
    spacing of 5 will result in a ball with radius 2 pixels.
    """
    radius_in_voxels = np.array([round(i) for i in radius / np.array(spacing)])
    n = 2 * radius_in_voxels + 1
    ball_iso = ball(max(n) * 2, dtype=np.float64)
    ball_resampled = resize(ball_iso, n, 1, 'constant', 0, clip=True, anti_aliasing=False, preserve_range=True)
    ball_resampled[ball_resampled > 0.5] = 1
    ball_resampled[ball_resampled <= 0.5] = 0
    return ball_resampled.astype(dtype)


def label_with_component_sizes(binary_image: np.ndarray, connectivity: int = None) -> Tuple[np.ndarray, dict]:
    if not binary_image.dtype == bool:
        print('Warning: it would be way faster if your binary image had dtype bool')
    labeled_image, num_components = label(binary_image, return_num=True, connectivity=connectivity)
    component_sizes = {i + 1: j for i, j in enumerate(np.bincount(labeled_image.ravel())[1:])}
    return labeled_image, component_sizes


def remove_all_but_largest_component(binary_image: np.ndarray, connectivity: int = None) -> np.ndarray:
    """
    Removes all but the largest component in binary_image. Replaces pixels that don't belong to it with background_label
    """
    filter_fn = lambda x, y: [i for i, j in zip(x, y) if j == max(y)]
    return generic_filter_components(binary_image, filter_fn, connectivity)


def generic_filter_components(binary_image: np.ndarray, filter_fn: Callable[[List[int], List[int]], List[int]],
                              connectivity: int = None):
    """
    filter_fn MUST return the component ids that should be KEPT!
    filter_fn will be called as: filter_fn(component_ids, component_sizes) and is expected to return a List of int

    returns a binary array that is True where the filtered components are
    """
    labeled_image, component_sizes = label_with_component_sizes(binary_image, connectivity)
    component_ids = list(component_sizes.keys())
    component_sizes = list(component_sizes.values())
    keep = filter_fn(component_ids, component_sizes)
    return np.in1d(labeled_image.ravel(), keep).reshape(labeled_image.shape)


def remove_components(binary_image: np.ndarray, threshold_size_in_pixels: int, threshold_type: str = 'min',
                      connectivity=None, verbose: bool = False):
    """
    This function removes either large or small components in a binary image, depending on the threshold_size_in_pixels.
    It is based on skimage's connected component labeling.
    If this function uses too much RAM, try out remove_components_cc3d.

    params:
        binary_image: 2D / 3D binary numpy image
        threshold_size_in_pixels: the threshold number of pixels/voxels of components to determine them as small / large
        threshold_type: 'min' means threshold_size_in_pixels defines minimal number of pixels -> removes small components
                        'max' means threshold_size_in_pixels defines maximal number of pixels -> removes large components
        connectivity: Maximum number of orthogonal hops to consider a pixel/voxel as a neighbor. Accepted values are
                      ranging from 1 to input.ndim. If None, a full connectivity of input.ndim is used.
        verbose: if True, prints number of removed components
    returns:
        binary image with removed components
    """

    assert threshold_type in ['min', 'max'], 'threshold type can only be min or max'
    binary_image = np.copy(binary_image)
    labeled_image, component_sizes = label_with_component_sizes(binary_image, connectivity)
    if threshold_type == 'min':
        keep = np.array([i for i, j in component_sizes.items() if j >= threshold_size_in_pixels])
        if verbose:
            print(f'{len(keep)} objects are larger than the minimum size of {threshold_size_in_pixels}. '
                  f'Removing {len(component_sizes) - len(keep)} small objects...')
    elif threshold_type == 'max':
        keep = np.array([i for i, j in component_sizes.items() if j < threshold_size_in_pixels])
        if verbose:
            print(f'{len(keep)} objects are smaller than the maximum size of {threshold_size_in_pixels}. '
                  f'Removing {len(component_sizes) - len(keep)} large objects...')

    keep = np.in1d(labeled_image, keep).astype(binary_image.dtype).reshape(binary_image.shape)
    return keep


def remove_components_cc3d(binary_image: np.ndarray, threshold_size_in_pixels: int, threshold_type: str = 'min',
                           connectivity=26, verbose: bool = False):
    """
    This function removes either large or small components in a binary image, depending on the threshold_size_in_pixels.
    It has a similar functionality as the skimage version remove_components but it uses connected-components-3d (https://github.com/seung-lab/connected-components-3d/) instead.
    In some use cases this function can be faster and use much less RAM than the skimage version.

    params:
        binary_image: 2D / 3D binary numpy image
        threshold_size_in_pixels: the threshold number of pixels/voxels of components to determine them as small / large
        threshold_type: 'min' means threshold_size_in_pixels defines minimal number of pixels -> removes small components
                        'max' means threshold_size_in_pixels defines maximal number of pixels -> removes large components
        connectivity: defines the neighborhood in which pixels count as connected, only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
                      2D:
                        - 4: neighbors that touch each others faces are counted (vertical and horizontal neighbors)
                        - 8: neighbors that touch each others faces or corners are counted (vertical, horizontal and diagonal neighbors)
                      3D:
                        - 6: neighbors that touch each others faces are counted
                        - 18: neighbors that touch each others faces or edges are counted
                        - 26: neighbors that touch each others faces, edges or corners are counted
        verbose: if True, prints number of removed components
    returns:
        binary image with removed components
    """

    assert threshold_type in ['min', 'max'], 'threshold type can only be min or max'

    if binary_image.ndim == 2:
        assert connectivity in [4, 8], 'Only connectivities of 4 or 8 are allowed for a 2D image'
    elif binary_image.ndim == 3:
        assert connectivity in [6, 18, 26], 'Only connectivities of 6, 18 or 26 are allowed for a 3D image'

    array = np.copy(binary_image)
    components, num_components = cc3d.connected_components(array, return_N=True, connectivity=connectivity)
    stats = cc3d.statistics(components)
    voxel_counts = stats['voxel_counts'][1:]  # first label is assumed to be background and is ignored

    if threshold_type == 'min':
        labels_to_keep = np.argwhere(voxel_counts >= threshold_size_in_pixels).reshape(-1) + 1
        if verbose:
            print(f'{len(labels_to_keep)} objects are larger than the minimum size of {threshold_size_in_pixels}. '
                  f'Removing {num_components - len(labels_to_keep)} small objects...')
    elif threshold_type == 'max':
        labels_to_keep = np.argwhere(voxel_counts < threshold_size_in_pixels).reshape(-1) + 1
        if verbose:
            print(f'{len(labels_to_keep)} objects are smaller than the maximum size of {threshold_size_in_pixels}. '
                  f'Removing {num_components - len(labels_to_keep)} large objects...')

    # in the component mask set every label that should be removed to 0
    keep = np.in1d(components, labels_to_keep).astype(binary_image.dtype).reshape(binary_image.shape)

    return keep


if __name__ == '__main__':
    print(generate_ball((10, 10, 10), (5, 5, 5)).shape)
    print(generate_ball((10, 5, 15), (1, 1, 1)).shape)
