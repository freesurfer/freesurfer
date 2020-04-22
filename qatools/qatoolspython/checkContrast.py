"""
This module provides a function to compute the WM/GM contrast SNR.

"""

# -----------------------------------------------------------------------------

def checkContrast(subjects_dir, subject):
    """
    A function to compute the WM/GM contrast SNR.

    This function computes the WM/GM contrast SNR based on the output of the 
    pctsurfcon function. 

    Required arguments: 
        - Subjects directory
        - Subject

    Returns:
        - con_lh_snr
        - con_rh_snr

    Requires surf/[lr]h.w-g.pct.mgh and label/[lr]h.cortex.label. If not found, 
    NaNs will be returned.

    """

    # Imports
    import os 
    import sys
    import numpy
    import nibabel
    from qatoolspython.qatoolspythonUtils import importMGH

    # Message
    print('Checking WM/GM contrast SNR ...')

    # Check if files exist
    path_pct_lh = os.path.join(subjects_dir,subject,"surf","lh.w-g.pct.mgh")
    if not os.path.exists(path_pct_lh):
        print("WARNING: could not find "+path_pct_lh+", returning NaNs")
        return numpy.nan

    path_pct_rh = os.path.join(subjects_dir,subject,"surf","rh.w-g.pct.mgh")
    if not os.path.exists(path_pct_rh):
        print("WARNING: could not find "+path_pct_rh+", returning NaNs")
        return numpy.nan

    path_label_cortex_lh  = os.path.join(subjects_dir,subject,"label","lh.cortex.label")
    if not os.path.exists(path_label_cortex_lh):
        print("WARNING: could not find "+path_label_cortex_lh+", returning NaNs")
        return numpy.nan

    path_label_cortex_rh = os.path.join(subjects_dir,subject,"label","rh.cortex.label")
    if not os.path.exists(path_label_cortex_rh):
        print("WARNING: could not find "+path_label_cortex_rh+", returning NaNs")
        return numpy.nan
    
    # Get the data fromt the mgh files
    con_lh = importMGH(path_pct_lh)
    con_rh = importMGH(path_pct_rh)
    
    label_array_lh = nibabel.freesurfer.io.read_label(path_label_cortex_lh)
    label_array_rh = nibabel.freesurfer.io.read_label(path_label_cortex_rh)    
    
    # Only take the values of the cortex to compute the contrast control
    con_lh = numpy.take(con_lh, label_array_lh)
    con_rh = numpy.take(con_rh, label_array_rh)
    
    # Compute the Contrast to noise ratio
    con_lh_mean = numpy.mean(con_lh)
    con_lh_std = numpy.std(con_lh)
    con_lh_snr = con_lh_mean/con_lh_std
    print("WM/GM contrast SNR for the left hemisphere:", '{:.4}'.format(con_lh_snr))
    
    con_rh_mean = numpy.mean(con_rh)
    con_rh_std = numpy.std(con_rh)
    con_rh_snr = con_rh_mean/con_rh_std 
    print("WM/GM contrast SNR for the right hemisphere:", '{:.4}'.format(con_rh_snr))
    
    # Return
    return con_lh_snr, con_rh_snr
    
