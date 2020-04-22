"""
This module provides a function to determine the rotation angles of the Talairach transform

"""

# -----------------------------------------------------------------------------

def checkRotation(subjects_dir, subject):
    """
    A function to determine the rotation angles of the Talairach transform.

    This function determines the rotation angles of the Talairach transform.
    It relies on the 'transforms3d' package by Matthew Brett 
    (https://matthew-brett.github.io/transforms3d).

    Required arguments:
        - subjects_dir : path to the subjects directory
        - subject : subject ID

    Returns:
        - rot_tal_x, rot_tal_y, rot_tal_z 

    Note: 
        - unit of rotations is radians

    Requires valid mri/transforms/talairach.lta file and the 'transforms3d' 
    python packge. If not found, NaNs will be returned.

    """

    #

    import os
    import numpy as np
    import importlib.util

    #

    if importlib.util.find_spec("transforms3d") is None:
        print('\nWARNING: the \'transforms3d\' package is required for running this script, returning NaNs.\n')
        return np.nan, np.nan, np.nan
    else:
        import transforms3d as tr

    # message

    print("Computing Talairach rotation angles ...")

    # read talairach.lta 

    if not os.path.isfile(os.path.join(subjects_dir,subject,"mri","transforms","talairach.lta")):
        print("WARNING: could not open "+os.path.join(subjects_dir,subject,"mri","transforms","talairach.lta")+", returning NaNs.")
        return np.nan, np.nan, np.nan

    with open(os.path.join(subjects_dir,subject,"mri","transforms","talairach.lta"), 'r') as datafile:
        lines = datafile.readlines()

    lines = lines[8:12]

    # convert to numpy array

    mat = list()
    for line in lines:
        mat.append([ float(x) for x in line.replace('\n','').replace(';','').split() ])
    mat = np.array(mat)

    # get translation, rotation, scale/zoom, and shear matrices

    T, R, Z, S = tr.affines.decompose44(mat)

    # now need to decompose R into euler angles; note this implies 'sxyz' 
    # (=static/extrinsic, 1:x, 2:y, 3:z rotation type and sequence)
    # https://matthew-brett.github.io/transforms3d/reference/transforms3d.euler.html

    rot_x, rot_y, rot_z = tr.euler.mat2euler(R)

    print('Found Talairach rotation angles: x =','{:.3}'.format(rot_x),', y =','{:.3}'.format(rot_y),', z =','{:.3}'.format(rot_z),'radians.')

    return rot_x, rot_y, rot_z
