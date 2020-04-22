"""
This module provides a function to check the relative size of the corpus callosum

"""

# -----------------------------------------------------------------------------

def checkCCSize(subjects_dir, subject):
    """
    A function to check the relative size of the corpus callosum.

    This function evaluates the relative size of the corpus callosum in order 
    to detect possible outliers. It computes the sum of the volumes of the 
    corpus callosum from in the aseg.stats file, and divides by total 
    intracranial volume.

    Required Arguments: 
        - subjects_dir
        - subject

    Returns:
        - relative_cc

    Requires a valid stats/aseg.stats file. If not found, NaNs are returned.

    """

    # Imports
    import os
    import numpy as np

    # Message
    print('Checking size of the corpus callosum ...')

    # Check if files exist
    path_stats_file = os.path.join(subjects_dir,subject,"stats","aseg.stats")

    try:
        with open(path_stats_file) as stats_file:
            aseg_stats = stats_file.read().splitlines()
    except FileNotFoundError:
        print("WARNING: could not open "+path_stats_file+", returning NaNs.")
        return np.nan

    # Initialize

    cc_elements = ['CC_Posterior', 'CC_Mid_Posterior', 'CC_Central', 'CC_Mid_Anterior', 'CC_Anterior']
    
    sum_cc = 0.0
    
    # Loop through the cc elements
    for cc_segmentation in cc_elements:
        # Look for the element in the aseg.stats line
        for aseg_stat_line in aseg_stats:
            # If the segmentation is found, compute the sum and return it. 
            if cc_segmentation in aseg_stat_line:
                sum_cc += float(aseg_stat_line.split()[3])
            elif 'EstimatedTotalIntraCranialVol' in aseg_stat_line:
                intracranial_volume= float(aseg_stat_line.split(',')[3])    

    relative_cc = sum_cc/intracranial_volume

    print("Relative size of the corpus callosum is",'{:.4}'.format(relative_cc))

    # Return
    return relative_cc
