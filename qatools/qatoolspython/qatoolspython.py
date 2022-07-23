"""
This module provides the main functionality of the qatoolspython package.

"""

# ==============================================================================
# FUNCTIONS

# ------------------------------------------------------------------------------
# get_help()

def get_help(print_help=True, return_help=False):
    """
    a function to return a help message

    """

    HELPTEXT="""

    qatools-python


    ============
    Description:
    ============

    This is a set of quality assurance / quality control scripts for Freesurfer 6.0
    processed structural MRI data.

    It is a revision, extension, and translation to the Python language of the 
    original Freesurfer QA Tools that are provided at 
    https://surfer.nmr.mgh.harvard.edu/fswiki/QATools

    It has been augmented by additional functions from the MRIQC toolbox, available 
    at https://github.com/poldracklab/mriqc and https://osf.io/haf97, and with code
    derived from the shapeDNA and brainPrint toolboxes, available at 
    https://reuter.mit.edu.

    The core functionality of this toolbox is to compute the following features:

    - wm_snr_orig   ...  signal-to-noise for white matter in orig.mgz
    - gm_snr_orig   ...  signal-to-noise for gray matter in orig.mgz
    - wm_snr_norm   ...  signal-to-noise for white matter in norm.mgz
    - gm_snr_norm   ...  signal-to-noise for gray matter in norm.mgz
    - cc_size       ...  relative size of the corpus callosum
    - holes_lh      ...  number of holes in the left hemisphere
    - holes_rh      ...  number of holes in the right hemisphere
    - defects_lh    ...  number of defects in the left hemisphere
    - defects_rh    ...  number of defects in the right hemisphere
    - topo_lh       ...  topological fixing time for the left hemisphere
    - topo_rh       ...  topological fixing time for the right hemisphere
    - con_snr_lh    ...  wm/gm contrast signal-to-noise ratio in the left hemisphere
    - con_snr_rh    ...  wm/gm contrast signal-to-noise ratio in the right hemisphere
    - rot_tal_x     ...  rotation component of the Talairach transform around the x axis
    - rot_tal_y     ...  rotation component of the Talairach transform around the y axis
    - rot_tal_z     ...  rotation component of the Talairach transform around the z axis

    The program will use an existing output directory (or try to create it) and 
    write a csv table into that location. The csv table will contain the above 
    metrics plus a subject identifier.
    
    In addition to the core functionality of the toolbox there are several optional
    modules that can be run according to need:

    - screenshots module

    This module allows for the automated generation of cross-sections of the brain 
    that are overlaid with the anatomical segmentations (asegs) and the white and 
    pial surfaces. These images will be saved to the 'screenshots' subdirectory 
    that will be created within the output directory. These images can be used for 
    quickly glimpsing through the processing results. Note that no display manager 
    is required for this module, i.e. it can be run on a remote server, for example.

    - fornix module

    This is a module to assess potential issues with the segmentation of the 
    corpus callosum, which may incorrectly include parts of the fornix. To assess 
    segmentation quality, a screesnhot of the contours of the corpus callosum 
    segmentation overlaid on the norm.mgz will be saved in the 'fornix' 
    subdirectory of the output directory. 
    
    - outlier module 
    
    This is a module to detect extreme values among the subcortical ('aseg') 
    segmentations. The outlier detection is based on comparisons with the 
    distributions of the sample as well as normative values taken from the 
    literature (see References).
    
    For comparisons with the sample distributions, extreme values are defined in 
    two ways: nonparametrically, i.e. values that are 1.5 times the interquartile 
    range below or above the 25th or 75th percentile of the sample, respectively, 
    and parametrically, i.e. values that are more than 2 standard deviations above 
    or below the sample mean. Note that a minimum of 5 supplied subjects is 
    required for running these analyses, otherwise `NaNs` will be returned.  
    
    For comparisons with the normative values, lower and upper bounds are computed 
    from the 95% prediction intervals of the regression models given in Potvin et 
    al., 1996, and values exceeding these bounds will be flagged. As an 
    alternative, users may specify their own normative values by using the 
    '--outlier-table' argument. This requires a custom csv table with headers
    `label`, `upper`, and `lower`, where `label` indicates a column of anatomical  
    names. It can be a subset and the order is arbitrary, but naming must exactly 
    match the nomenclature of the 'aseg.stats' file. `upper` and `lower` are user-
    specified upper and lower bounds.    
      
    The main csv table will be appended with the following summary variables, and 
    more detailed output about will be saved as csv tables in the 'outliers' 
    subdirectory of the main output directory.
    
    n_outliers_sample_nonpar ... number of structures that are 1.5 times the IQR  
                                 above/below the 75th/25th percentile  
    n_outliers_sample_param  ... number of structures that are 2 SD above/below 
                                 the mean
    n_outliers_norms         ... number of structures exceeding the upper and      
                                 lower bounds of the normative values
    

    ======
    Usage: 
    ======

        python3 qatools.py --subjects_dir <directory> --output_dir <directory>
                                  [--subjects SubjectID [SubjectID ...]] 
                                  [--screenshots] [--fornix] [-h]

        required arguments:
          --subjects_dir <directory>
                                subjects directory with a set of Freesurfer 6.0 
                                processed individual datasets.
          --output_dir <directory>
                                output directory

        optional arguments:
          --subjects SubjectID [SubjectID ...]
                                list of subject IDs
          --screenshots         create screenshots of individual brains
          --fornix              check fornix segmentation
          --outlier             run outlier detection
          --outlier-table       specify normative values (only in conjunction with
                                --outlier)          

        getting help:
          -h, --help            display this help message and exit


    ========================
    Use as a python package:
    ========================

    As an alternative to their command-line usage, the qc scripts can also be run 
    within a pure python environment, i.e. installed and imported as a python 
    package. 

    Use `import qatoolspython` (or sth. equivalent) to import the package within a 
    python environment.
    
    Use the `run_qatools` function from the `qatoolspython` module to run an analysis:
    
    `from qatoolspython import qatoolspython`
    
    `qatoolspython.run_qatools(subjects_dir='/my/subjects/dir',output_dir='/my/output/dir')`
    
    See `help(qatoolspython)` for further usage info and options.


    =============
    Known Issues: 
    =============

    The program will analyze recon-all logfiles, and may fail or return erroneous
    results if the logfile is append by multiple restarts of recon-all runs. 
    Ideally, the logfile should therefore consist of just a single, successful 
    recon-all run.


    ========
    Authors: 
    ========

    - qatools-python: Kersten Diers, Tobias Wolff, and Martin Reuter.
    - Freesurfer QA Tools: David Koh, Stephanie Lee, Jenni Pacheco, Vasanth Pappu, 
      and Louis Vinke. 
    - shapeDNA and brainPrint toolboxes: Martin Reuter


    ===========
    References:
    ===========

    Esteban O, Birman D, Schaer M, Koyejo OO, Poldrack RA, Gorgolewski KJ; MRIQC: 
    Advancing the Automatic Prediction of Image Quality in MRI from Unseen Sites; 
    PLOS ONE 12(9):e0184661; doi:10.1371/journal.pone.0184661.

    Wachinger C, Golland P, Kremen W, Fischl B, Reuter M; 2015; BrainPrint: a 
    Discriminative Characterization of Brain Morphology; Neuroimage: 109, 232-248; 
    doi:10.1016/j.neuroimage.2015.01.032.

    Reuter M, Wolter FE, Shenton M, Niethammer M; 2009; Laplace-Beltrami 
    Eigenvalues and Topological Features of Eigenfunctions for Statistical Shape 
    Analysis; Computer-Aided Design: 41, 739-755, doi:10.1016/j.cad.2009.02.007.

    Potvin O, Mouiha A, Dieumegarde L, Duchesne S, & Alzheimer's Disease Neuroimaging 
    Initiative; 2016; Normative data for subcortical regional volumes over the lifetime 
    of the adult human brain; Neuroimage: 137, 9-20; 
    doi.org/10.1016/j.neuroimage.2016.05.016

    =============
    Requirements:
    =============

    A working installation of Freesurfer 6.0 or later must be sourced.

    At least one subject whose structural MR image was processed with Freesurfer 
    6.0 or later.

    A Python version >= 3.5 is required to run this script.

    Required packages include (among others) the nibabel and skimage package for 
    the core functionality, plus the the matplotlib, pandas, and transform3d 
    packages for some optional functions and modules. See the `requirements.txt` 
    file for a complete list. Use `pip3 install -r requirements.txt` to install 
    these packages.
    
    This software has been tested on Ubuntu 16.04, CentOS7, and MacOS 10.14. 


    ========
    License:
    ========

    This software is licensed under the MIT License, see associated LICENSE file 
    for details.

    Copyright (c) 2019 Image Analysis Group, DZNE e.V.

    """

    if print_help:
        print(HELPTEXT)

    if return_help:
        return HELPTEXT


# ------------------------------------------------------------------------------
# parse_arguments

def _parse_arguments():
    """
    an internal function to parse input arguments

    """

    # imports
    import sys
    import argparse

    # parse
    parser = argparse.ArgumentParser(description='''
        This program takes existing Freesurfer 6.0 analysis results of one
        or more subjects and computes a set of quality metrics. These will be 
        reported in a summary csv table.

        For a description of these metrics, see the gitlab/github page or the 
        header section of this script.

        The (optional) analysis of shape features requires additional scripts 
        that can be obtained from https://reuter.mit.edu
        ''', 
        add_help=False, formatter_class=argparse.RawTextHelpFormatter)

    required = parser.add_argument_group('required arguments')
    required.add_argument('--subjects_dir', dest="subjects_dir", help="subjects directory with a set of Freesurfer 6.0 \nprocessed individual datasets.", metavar="<directory>", required=True)
    required.add_argument('--output_dir', dest="output_dir", help="output directory", metavar="<directory>", required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--subjects', dest="subjects", help="list of subject IDs. If omitted, all suitable sub-\ndirectories witin the subjects directory will be \nused.", default=[], nargs='+', metavar="SubjectID", required=False)
    #optional.add_argument('--shape', dest='shape', help="run shape analysis (requires additional scripts)", default=False, action="store_true", required=False)
    optional.add_argument('--shape', dest='shape', help=argparse.SUPPRESS, default=False, action="store_true", required=False) # shape is currently a hidden option
    optional.add_argument('--screenshots', dest='screenshots', help="create screenshots of individual brains", default=False, action="store_true", required=False)
    optional.add_argument('--fornix', dest='fornix', help="check fornix segmentation", default=False, action="store_true", required=False)
    optional.add_argument('--outlier', dest='outlier', help="run outlier detection", default=False, action="store_true", required=False)
    optional.add_argument('--outlier-table', dest="outlier_table", help="specify normative values", default=None, metavar="<filename>", required=False)

    help = parser.add_argument_group('getting help')
    help.add_argument('-h', '--help', help="display this help message and exit", action='help')

    # check if there are any inputs; if not, print help and exit
    if len(sys.argv)==1:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    return args.subjects_dir, args.output_dir, args.subjects, args.shape, args.screenshots, args.fornix, args.outlier, args.outlier_table

# ------------------------------------------------------------------------------
# check arguments

def _check_arguments(subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table):
    """
    an internal function to check input arguments

    """

    # --------------------------------------------------------------------------
    # imports

    import os
    import sys
    import errno

    import tempfile
    import importlib.util

    # --------------------------------------------------------------------------
    # check arguments

    # check if subject directory exists
    if os.path.isdir(subjects_dir):
        print("Found subjects directory", subjects_dir)
    else:
        print('ERROR: subjects directory '+subjects_dir+' is not an existing directory\n')
        sys.exit(1)

    # check if output directory exists or can be created and is writable
    if os.path.isdir(output_dir):
        print("Found output directory", output_dir)
    else:
        try:
            os.mkdir(output_dir)
        except:
            print('ERROR: cannot create output directory '+output_dir+'\n')
            sys.exit(1)

        try:
            testfile = tempfile.TemporaryFile(dir=output_dir)
            testfile.close()
        except OSError as e:
            if e.errno != errno.EACCES:  # 13
                e.filename = output_dir
                raise
            print('\nERROR: '+output_dir+' not writeable (check access)!\n')
            sys.exit(1)

    # check if screenshots subdirectory exists or can be created and is writable
    if screenshots is True:
        if os.path.isdir(os.path.join(output_dir,'screenshots')):
            print("Found screenshots directory", os.path.join(output_dir,'screenshots'))
        else:
            try:
                os.mkdir(os.path.join(output_dir,'screenshots'))
            except:
                print('ERROR: cannot create screenshots directory '+os.path.join(output_dir,'screenshots')+'\n')
                sys.exit(1)

            try:
                testfile = tempfile.TemporaryFile(dir=os.path.join(output_dir,'screenshots'))
                testfile.close()
            except OSError as e:
                if e.errno != errno.EACCES:  # 13
                    e.filename = os.path.join(output_dir,'screenshots')
                    raise
                print('\nERROR: '+os.path.join(output_dir,'screenshots')+' not writeable (check access)!\n')
                sys.exit(1)

    # check further screenshots dependencies
    if screenshots is True and importlib.util.find_spec("pandas") is None:
        print('\nERROR: the \'pandas\' package is required for running this script, please install.\n')
        sys.exit(1)

    if screenshots is True and importlib.util.find_spec("matplotlib") is None:
        print('\nERROR: the \'matplotlib\' package is required for running this script, please install.\n')
        sys.exit(1)

    if screenshots is True:
        path_check = os.path.join(os.environ.get('FREESURFER_HOME'), 'FreeSurferColorLUT.txt')
        if not os.path.isfile(path_check):
            print('\nERROR: the \'FreeSurferColorLUT.txt\' file needs to be present in the FREESURFER_HOME directory.\n')
            sys.exit(1)

    # check if fornix subdirectory exists or can be created and is writable
    if fornix is True:
        if os.path.isdir(os.path.join(output_dir,'fornix')):
            print("Found fornix directory", os.path.join(output_dir,'fornix'))
        else:
            try:
                os.mkdir(os.path.join(output_dir,'fornix'))
            except:
                print('ERROR: cannot create fornix directory '+os.path.join(output_dir,'fornix')+'\n')
                sys.exit(1)

            try:
                testfile = tempfile.TemporaryFile(dir=os.path.join(output_dir,'fornix'))
                testfile.close()
            except OSError as e:
                if e.errno != errno.EACCES:  # 13
                    e.filename = os.path.join(output_dir,'fornix')
                    raise
                print('\nERROR: '+os.path.join(output_dir,'fornix')+' not writeable (check access)!\n')
                sys.exit(1)

    # check if shape subdirectory exists or can be created and is writable
    if shape is True:
        if os.path.isdir(os.path.join(output_dir, 'brainprint')):
            print("Found brainprint directory", os.path.join(output_dir,'brainprint'))
        else:
            try:
                os.makedirs(os.path.join(output_dir,'brainprint'))
            except:
                print('\nERROR: cannot create brainprint directory '+os.path.join(output_dir, 'brainprint')+'\n')
                sys.exit(1)

            try:
                testfile = tempfile.TemporaryFile(dir=os.path.join(output_dir,'brainprint'))
                testfile.close()
            except OSError as e:
                if e.errno != errno.EACCES:  # 13
                    e.filename = os.path.join(output_dir,'brainprint')
                    raise
                print('\nERROR: '+os.path.join(output_dir,'brainprint')+' not writeable (check access)!\n')
                sys.exit(1)

    # check if shapeDNA / brainPrint dependencies
    if shape is True:
        # check if brainprintpython can be imported
        if  importlib.util.find_spec("brainprintpython") is None:
            print("\nERROR: could not import the brainprintpython package, is it installed?") 
            sys.exit(1)

    # check if outlier subdirectory exists or can be created and is writable
    if outlier is True:
        if os.path.isdir(os.path.join(output_dir, 'outliers')):
            print("Found outliers directory", os.path.join(output_dir,'outliers'))
        else:
            try:
                os.makedirs(os.path.join(output_dir,'outliers'))
            except:
                print('\nERROR: cannot create outliers directory '+os.path.join(output_dir, 'outliers')+'\n')
                sys.exit(1)

            try:
                testfile = tempfile.TemporaryFile(dir=os.path.join(output_dir,'outliers'))
                testfile.close()
            except OSError as e:
                if e.errno != errno.EACCES:  # 13
                    e.filename = os.path.join(output_dir,'outliers')
                    raise
                print('\nERROR: '+os.path.join(output_dir,'outliers')+' not writeable (check access)!\n')
                sys.exit(1)

    # check if outlier-table exists if it was given, otherwise exit
    if outlier_table is not None:
        if os.path.isfile(outlier_table):
            print("Found table with normative values ", outlier_table)
        else:
            print("Could not find table with normative values ", outlier_table)
            sys.exit(1)

    # if subjects are not given, get contents of the subject directory and 
    # check if aseg.stats (as a proxy) exists
    if subjects == []:
        for subject in os.listdir(subjects_dir):
            path_aseg_stat = os.path.join(subjects_dir,subject,"stats","aseg.stats")
            if os.path.isfile(path_aseg_stat):
                print("Found subject",subject)
                subjects.extend([subject])

    # check for required files
    subjects_to_remove = list()
    for subject in subjects:

        # -files: stats/aseg.stats
        path_check = os.path.join(subjects_dir,subject,"stats","aseg.stats")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        # -files: surf/[lr]h.w-g.pct.mgh, label/[lr]h.cortex.label
        path_check = os.path.join(subjects_dir, subject, "surf", "lh.w-g.pct.mgh")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        path_check = os.path.join(subjects_dir, subject, "surf", "rh.w-g.pct.mgh")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        path_check = os.path.join(subjects_dir, subject, "label", "lh.cortex.label")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        path_check = os.path.join(subjects_dir, subject, "label", "rh.cortex.label")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        # -files: mri/transforms/talairach.lta
        path_check = os.path.join(subjects_dir, subject, "mri", "transforms", "talairach.lta")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        # -files: mri/norm.mgz, mri/aseg.mgz, mri/aparc+aseg.mgz
        path_check = os.path.join(subjects_dir, subject, "mri", "norm.mgz")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        path_check = os.path.join(subjects_dir, subject, "mri", "aseg.mgz")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        path_check = os.path.join(subjects_dir, subject, "mri", "aparc+aseg.mgz")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        # -files: scripts/recon-all.log
        path_check = os.path.join(subjects_dir, subject, "scripts", "recon-all.log")
        if not os.path.isfile(path_check):
            print("Could not find", path_check, "for subject", subject)
            subjects_to_remove.extend([subject])

        # check screenshots
        if screenshots is True:

            # -files: surf/[lr]h.white (optional), surf/[lr]h.pial (optional)
            path_check = os.path.join(subjects_dir, subject, "surf", "lh.white")
            if not os.path.isfile(path_check):
                print("Could not find", path_check, "for subject", subject)
                subjects_to_remove.extend([subject])

            path_check = os.path.join(subjects_dir, subject, "surf", "rh.white")
            if not os.path.isfile(path_check):
                print("Could not find", path_check, "for subject", subject)
                subjects_to_remove.extend([subject])

            path_check = os.path.join(subjects_dir, subject, "surf", "lh.pial")
            if not os.path.isfile(path_check):
                print("Could not find", path_check, "for subject", subject)
                subjects_to_remove.extend([subject])

            path_check = os.path.join(subjects_dir, subject, "surf", "rh.pial")
            if not os.path.isfile(path_check):
                print("Could not find", path_check, "for subject", subject)
                subjects_to_remove.extend([subject])

        # check fornix
        if fornix is True:

            # -files: mri/transforms/cc_up.lta
            path_check = os.path.join(subjects_dir, subject, "mri", "transforms", "cc_up.lta")
            if not os.path.isfile(path_check):
                print("Could not find", path_check, "for subject", subject)
                subjects_to_remove.extend([subject])

    # remove subjects with missing files
    [ subjects.remove(x) for x in subjects_to_remove ]

    # check if we have any subjects after all
    if subjects == []:
        print("\nERROR: no subjects to process") 
        sys.exit(1)

    # now return
    return subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table


# ------------------------------------------------------------------------------
# check packages

def _check_packages():
    """
    an internal function to check required / recommended packages

    """

    import os
    import sys
    import importlib.util

    if os.environ.get('FREESURFER_HOME') is None:
        print('\nERROR: need to set the FREESURFER_HOME environment variable\n')
        sys.exit(1)

    if sys.version_info <= (3, 5):
        print('\nERROR: Python version must be 3.5 or greater\n')
        sys.exit(1)

    if importlib.util.find_spec("skimage") is None:
        print('\nERROR: the \'skimage\' package is required for running this script, please install.\n')
        sys.exit(1)

    if importlib.util.find_spec("nibabel") is None:
        print('\nERROR: the \'nibabel\' package is required for running this script, please install.\n')
        sys.exit(1)

    if importlib.util.find_spec("transforms3d") is None:
        # this package is less important and less standard, so we just return a
        # warning (and NaNs) if it is not found.
        print('\nWARNING: the \'transforms3d\' package is recommended, please install.\n')


# ------------------------------------------------------------------------------
# do qatools

def _do_qatools(subjects_dir, output_dir, subjects, shape=False, screenshots=False, fornix=False, outlier=False, outlier_table=None):
    """
    an internal function to run the qatools submodules

    """

    # ------------------------------------------------------------------------------
    # imports

    import os
    import csv
    import time

    from qatoolspython.checkSNR import checkSNR
    from qatoolspython.checkCCSize import checkCCSize
    from qatoolspython.checkTopology import checkTopology
    from qatoolspython.checkContrast import checkContrast
    from qatoolspython.checkRotation import checkRotation
    from qatoolspython.evaluateFornixSegmentation import evaluateFornixSegmentation
    from qatoolspython.createScreenshots import createScreenshots
    from qatoolspython.outlierDetection import outlierTable
    from qatoolspython.outlierDetection import outlierDetection

    # ------------------------------------------------------------------------------
    # internal settings (might be turned into command-line arguments in the future)

    SNR_AMOUT_EROSION = 3
    FORNIX_SCREENSHOT = True
    FORNIX_SHAPE = False
    FORNIX_N_EIGEN = 15
    OUTLIER_N_MIN = 5

    # --------------------------------------------------------------------------
    # process

    # start the processing with a message
    print("")
    print("-----------------------------")

    # create dict for this subject
    metricsDict = dict()

    # loop through the specified subjects
    for subject in subjects:

        #
        print("Starting qatools-python for subject", subject, "at", time.strftime('%Y-%m-%d %H:%M %Z', time.localtime(time.time())))
        print("")

        # ----------------------------------------------------------------------
        # compute core metrics

        # get WM and GM SNR for orig.mgz
        wm_snr_orig, gm_snr_orig = checkSNR(subjects_dir, subject, SNR_AMOUT_EROSION, ref_image="orig.mgz")

        # get WM and GM SNR for norm.mgz
        wm_snr_norm, gm_snr_norm = checkSNR(subjects_dir, subject, SNR_AMOUT_EROSION, ref_image="norm.mgz")

        # check CC size
        cc_size = checkCCSize(subjects_dir, subject)

        # check topology
        holes_lh, holes_rh, defects_lh, defects_rh, topo_lh, topo_rh = checkTopology(subjects_dir, subject)

        # check contrast
        con_snr_lh, con_snr_rh = checkContrast(subjects_dir, subject)

        # check rotation
        rot_tal_x, rot_tal_y, rot_tal_z = checkRotation(subjects_dir, subject)

        # store data
        metricsDict.update( { subject : {
            'subject' : subject,
            'wm_snr_orig': wm_snr_orig, 'gm_snr_orig' : gm_snr_orig, 
            'wm_snr_norm' : wm_snr_norm, 'gm_snr_norm' : gm_snr_norm, 
            'cc_size' : cc_size, 
            'holes_lh' : holes_lh, 'holes_rh' : holes_rh, 'defects_lh' : defects_lh, 'defects_rh' : defects_rh, 'topo_lh' : topo_lh, 'topo_rh' : topo_rh,
            'con_snr_lh' : con_snr_lh, 'con_snr_rh' : con_snr_rh,
            'rot_tal_x' : rot_tal_x, 'rot_tal_y' : rot_tal_y , 'rot_tal_z' : rot_tal_z  
            }})

        # 
        print("")

        # ----------------------------------------------------------------------
        # run optional modules: shape analysis

        if shape is True:
    
            # message
            print("-----------------------------")
            print("Running brainPrint analysis ...")
            print("")

            # compute brainprint (will also compute shapeDNA)
            from brainprintpython import pyBrainPrint
            from brainprintpython import pyPostProc

            # check / create subject-specific brainprint_outdir
            brainprint_outdir = os.path.join(output_dir,'brainprint',subject)
            if not os.path.isdir(brainprint_outdir):
                os.mkdir(brainprint_outdir)

            # set options
            class options:
                sdir = subjects_dir
                sid = subject
                outdir = brainprint_outdir
                evec = True
                skipcortex = True
                num = 15
                bcond = 1
                brainprint = os.path.join(brainprint_outdir, subject+'.brainprint.csv')

            class optionsPostProc:
                file = os.path.join(brainprint_outdir, subject+'.brainprint.csv')
                csvfiles = [ os.path.join(brainprint_outdir, subject+'.brainprint.csv') ]
                out =  brainprint_outdir
                outcov = None
                vol = 1
                lin = True
                asy = "euc"

            # run brainPrint
            structures, evmat = pyBrainPrint.compute_brainprint(options)

            # write EVs
            pyBrainPrint.write_ev(options, structures, evmat)

            # run postProc
            postProcDict = pyPostProc.compute_postproc(optionsPostProc)

            # get a subset of the brainprint results
            distDict = { subject : postProcDict[os.path.join(brainprint_outdir, subject+".brainprint.csv")]['dist'] }
    
            # store data
            metricsDict[subject].update(distDict[subject])

        # ----------------------------------------------------------------------
        # run optional modules: screenshots

        if screenshots is True:

            # message
            print("-----------------------------")
            print("Creating screenshots ...")
            print("")

            # check / create subject-specific screenshots_outdir
            screenshots_outdir = os.path.join(output_dir,'screenshots',subject)
            if not os.path.isdir(screenshots_outdir):
                os.mkdir(screenshots_outdir)
            outfile = os.path.join(screenshots_outdir,subject+'.png')

            # process
            createScreenshots(SUBJECT=subject, SUBJECTS_DIR=subjects_dir, OUTFILE=outfile, INTERACTIVE=False)

        # ----------------------------------------------------------------------
        # run optional modules: fornix

        if fornix is True:

            # message
            print("-----------------------------")
            print("Checking fornix segmentation ...")
            print("")

            # check / create subject-specific fornix_outdir
            fornix_outdir = os.path.join(output_dir,'fornix',subject)
            if not os.path.isdir(fornix_outdir):
                os.mkdir(fornix_outdir)

            # process
            fornixShapeOutput = evaluateFornixSegmentation(SUBJECT=subject,SUBJECTS_DIR=subjects_dir,OUTPUT_DIR=fornix_outdir,CREATE_SCREENSHOT=FORNIX_SCREENSHOT,RUN_SHAPEDNA=FORNIX_SHAPE,N_EIGEN=FORNIX_N_EIGEN)

            # create a dictionary from fornix shape ouput
            fornixShapeDict = { subject : dict(zip(map("fornixShapeEV{:0>3}".format,range(FORNIX_N_EIGEN)), fornixShapeOutput)) }
       
            # store data
            if FORNIX_SHAPE:
                metricsDict[subject].update(fornixShapeDict[subject])

        # message
        print("Finished subject", subject, "at", time.strftime('%Y-%m-%d %H:%M %Z', time.localtime(time.time())))
        print("")

    # --------------------------------------------------------------------------
    # run optional modules: outlier detection

    if outlier is True:

        # message
        print("---------------------------------------")
        print("Running outlier detection module ...")
        print("")

        # determine outlier-table and get data
        if outlier_table is None:
            outlierDict = outlierTable()
        else:
            outlierDict = dict()
            with open(outlier_table, newline='') as csvfile:
                outlierCsv = csv.DictReader(csvfile, delimiter=',')
                for row in outlierCsv:
                    outlierDict.update({row['label']: {'lower': float(row['lower']), 'upper': float(row['upper'])}})

        # process
        outlier_outdir = os.path.join(output_dir, 'outliers')
        n_outlier_sample_nonpar, n_outlier_sample_param, n_outlier_norms = outlierDetection(subjects, subjects_dir, outlier_outdir, outlierDict, min_no_subjects=OUTLIER_N_MIN)

        # create a dictionary from outlier module ouput
        outlierDict = dict()
        for subject in subjects:
            outlierDict.update({subject : {
                'n_outlier_sample_nonpar' : n_outlier_sample_nonpar[subject],
                'n_outlier_sample_param': n_outlier_sample_param[subject],
                'n_outlier_norms': n_outlier_norms[subject]
                }
            })

        # store data
        for subject in subjects:
            metricsDict[subject].update(outlierDict[subject])

        # message
        print("Done")
        print("")

    # --------------------------------------------------------------------------
    # generate output

    # we pre-specify the fieldnames because we want to have this particular order
    metricsFieldnames = ['subject','wm_snr_orig','gm_snr_orig','wm_snr_norm','gm_snr_norm','cc_size','holes_lh','holes_rh','defects_lh','defects_rh','topo_lh','topo_rh','con_snr_lh','con_snr_rh','rot_tal_x', 'rot_tal_y', 'rot_tal_z']

    if shape is True:
        metricsFieldnames.extend(distDict[subject].keys())

    if fornix is True and FORNIX_SHAPE is True:
        metricsFieldnames.extend(sorted(fornixShapeDict[subject].keys()))

    if outlier is True:
        metricsFieldnames.extend(sorted(outlierDict[subject].keys()))

    # determine output file names
    path_data_file = os.path.join(output_dir,'qatools-results.csv')

    # write csv
    with open(path_data_file, 'w') as datafile:
        csvwriter = csv.DictWriter(datafile, fieldnames=metricsFieldnames, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writeheader()
        for subject in sorted(list(metricsDict.keys())):
            csvwriter.writerow(metricsDict[subject])


# ------------------------------------------------------------------------------
# run qatools

def run_qatools(subjects_dir, output_dir, subjects=[], shape=False, screenshots=False, fornix=False, outlier=False, outlier_table=None):
    """
    a function to run the qatools submodules

    """

    # ------------------------------------------------------------------------------
    #

    # check arguments
    subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table = _check_arguments(subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table)

    # check packages
    _check_packages()

    # run qatools
    _do_qatools(subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table)
