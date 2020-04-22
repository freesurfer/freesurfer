"""
This module provides a function to check the topology of left and right surfaces

"""

# -----------------------------------------------------------------------------

def checkTopology(subjects_dir, subject):
    """
    A function to check the topology of left and right surfaces.

    This scripts extract the information about the number of holes and defects 
    in the left and right hemisphere, and also returns topological fixing time.

    Required arguments:
        - Subjects directory
        - Subject

    Returns:
        - lh_holes, rh_holes, lh_defects, rh_defects, topo_time_lh, topo_time_rh

    Requires valid scripts/recon-all.log file. If not found, NaNs will be 
    returned.

    """

    # Imports

    import os
    import numpy as np

    # Message

    print("Checking topology of the surfaces ...")

    # Get the logfile, and return with NaNs if unsuccessful:

    path_log_file = os.path.join(subjects_dir,subject,"scripts","recon-all.log")

    try:
        with open(path_log_file, 'r') as logfile:
            lines_log_file = logfile.read().splitlines()
    except FileNotFoundError:
            print("WARNING: could not find "+path_log_file+", returning NaNs.")
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    # Extract info from logfile

    foundDefectsLH = False
    foundTopoLH = False

    for line_log_file in lines_log_file:

        # Look for the number of holes in the left and right hemisphere
        if "orig.nofix lhholes" in line_log_file:
            lh_holes = line_log_file.split()[3]
            lh_holes = lh_holes[:-1]
            lh_holes = int(lh_holes)
            rh_holes = line_log_file.split()[6]
            print("Number of holes in the left hemisphere:", lh_holes)
            print("Number of holes in the right hemisphere:", rh_holes)

        # Look for the number of defects
        if "defects found" in line_log_file and foundDefectsLH is False :
            lh_defects = line_log_file.split()[0]
            lh_defects = int(lh_defects)
            print("Number of defects in the left hemisphere:", lh_defects)
            foundDefectsLH = True
        elif "defects found" in line_log_file and foundDefectsLH is True :
            rh_defects = line_log_file.split()[0]
            rh_defects = int(rh_defects)
            print("Number of defects in the right hemisphere:", rh_defects)

        # Look for the topological fixing time in the log file
        if "topology fixing took" in line_log_file and foundTopoLH is False:
            topo_time_lh = line_log_file.split()[3]
            print("Topological fixing time for the left hemisphere:", topo_time_lh, "min")
            foundTopoLH = True
        elif "topology fixing took" in line_log_file and foundTopoLH is True:
            topo_time_rh = line_log_file.split()[3]
            print("Topological fixing time for the right hemisphere:", topo_time_rh, "min")

    # Return

    return lh_holes, rh_holes, lh_defects, rh_defects, topo_time_lh, topo_time_rh
