import numpy as np

from .deprecations import deprecate, replace, unsure, notneeded


@deprecate('Im guessing no one uses this directly, but if you do check out sf.transform.orientation')
def build_matrix(matrix, shape=(256, 256, 256), voxsize=(1, 1, 1)):
    matrix[:3, :3] *= voxsize  # TODO this needs to be copied
    if shape is not None:
        pcrs = np.append((np.array(shape) - 1) / 2, 1)
        matrix[:3, 3] = -np.matmul(matrix, pcrs)[:3]
    return matrix


@deprecate('Im guessing no one uses this directly, but if you do check out sf.transform.orientation')
def orientation_from_matrix(matrix):
    """
    Examines a 4x4 direction cosine matrix and determines the corresponding orientation
    string, which indicates the primary direction of each axis in the matrix. The string
    characters can be L/R, A/P, I/S, which correspond to Left/Right, Anterior/Posterior,
    and Inferior/Superior, respectively.
    """
    matrix = matrix[:3, :3]
    orientation = ''
    for i in range(3):
        sag, cor, ax = matrix[:, i]
        if np.abs(sag) > np.abs(cor) and np.abs(sag) > np.abs(ax):
            orientation += 'R' if sag > 0 else 'L'
        elif np.abs(cor) > np.abs(ax):
            orientation += 'A' if cor > 0 else 'P'
        else:
            orientation += 'S' if ax > 0 else 'I'
    return orientation


@deprecate('Im guessing no one uses this directly, but if you do check out sf.transform.orientation')
def matrix_from_orientation(orientation):
    """
    Computes the 4x4 direction cosine matrix corresponding to an orientation string.
    """
    orientation = orientation.upper()
    check_orientation(orientation)

    matrix = np.zeros((4, 4))
    for i, c in enumerate(orientation):
        matrix[:3, i] -= [c == x for x in 'LPI']
        matrix[:3, i] += [c == x for x in 'RAS']
    matrix[3, 3] = 1
    return matrix


@deprecate('Im guessing no one uses this directly, but if you do check out sf.transform.orientation')
def check_orientation(orientation):
    """
    Checks an orientation string to ensure it is valid, meaning that all axes are represented
    exactly once and no invalid characters are present. Throws `ValueError` if orientation
    is invalid.
    """
    orientation = orientation.upper()

    if len(orientation) != 3:
        raise ValueError('Bad orientation: expected 3 characters, got %d.' % len(orientation))

    axes = ['LR', 'PA', 'IS']
    rs = np.zeros(3, dtype='int')
    for c in orientation:
        r = [c in axis for axis in axes]
        if not any(r):
            raise ValueError('Bad orientation: unknown character "%s".' % c)
        rs += r

    for i in range(3):
        if rs[i] > 1:
            raise ValueError('Bad orientation: %s axis represented multiple times.' % axes[i])
        if rs[i] == 0:
            raise ValueError('Bad orientation: %s axis not represented.' % axes[i])


@deprecate('Im guessing no one uses this directly, but if you do check out sf.transform.orientation')
def slice_direction(orientation):
    """
    Determines the primary slice plane from an orientation string.
    """
    check_orientation(orientation)
    if orientation[2] in 'LR':
        return 'sagittal'
    elif orientation[2] in 'PA':
        return 'coronal'
    elif orientation[2] in 'IS':
        return 'axial'
    else:
        raise ValueError('Invalid orientation: %s.' % orientation)
