"""
This module provides a function to evaluate potential missegmentation of the fornix

"""

# -----------------------------------------------------------------------------

def evaluateFornixSegmentation(SUBJECT, SUBJECTS_DIR, OUTPUT_DIR, CREATE_SCREENSHOT = True, RUN_SHAPEDNA = True, N_EIGEN = 15):
    """
    A function to evaluate potential missegmentation of the fornix.

    This script evaluates potential missegmentation of the fornix, which may 
    erroneously be attached to the 'corpus collosum' label.

    It will run Freesurfer's 'mri_convert' script to apply the cc_up.lta 
    transform to the norm.mgz and the aseg files, and create a binary corpus 
    callosum mask and surface. Resulting files are saved to subject-specific
    directories witin the 'fornix' subdirectory of the output directory.

    If the corresponding arguments are set to 'True', the script will also 
    create screenshots and run a shape analysis of the corpus callosum surface.
    Resulting files will be saved to the same directory as indicated above.

    Required arguments:
        - SUBJECT
        - SUBJECTS_DIR
        - OUTPUT_DIR

    Optional arguments:
        - CREATE_SCREENSHOT <bool> (default: True)
        - RUN_SHAPEDNA <bool> (default: True)
        - N_EIGEN <int> number of Eigenvalues for shape analyis (default: 30)

    Requires (if not found, returns NaNs):
        - mri/transforms/cc_up.lta
        - mri/aseg.mgz
        - mri/norm.mgz

    Returns:
        - a numpy array of N_EIGEN eigenvalues if RUN_SHAPEDNA is True, 
          otherwise a numpy array of NaNs of the same dimension

    """

    # --------------------------------------------------------------------------
    # imports

    import os
    import sys
    import shlex
    import subprocess
    import numpy as np
    from qatoolspython.createScreenshots import createScreenshots

    # --------------------------------------------------------------------------
    # auxiliary functions

    def split_callback(option, opt, value, parser):
      setattr(parser.values, option.dest, value.split(','))

    def which(program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
            if is_exe(os.path.join(os.getenv('SHAPEDNA_HOME'),program)):
                return os.path.join(os.getenv('SHAPEDNA_HOME'),program)
            if is_exe(os.path.join('.',program)):
                return os.path.join('.',program)

        return None

    def my_print(message):
        """
        print message, then flush stdout
        """
        print(message)
        sys.stdout.flush()

    def run_cmd(cmd,err_msg):
        """
        execute the comand
        """
        clist = cmd.split()
        progname=which(clist[0])
        if (progname) is None:
            my_print('ERROR: '+ clist[0] +' not found in path!')
            sys.exit(1)
        clist[0]=progname
        cmd = ' '.join(clist)
        my_print('#@# Command: ' + cmd+'\n')

        args = shlex.split(cmd)
        try:
            subprocess.check_call(args)
        except subprocess.CalledProcessError as e:
            my_print('ERROR: '+err_msg)
            #sys.exit(1)
            raise
        my_print('\n')

    # --------------------------------------------------------------------------
    # main part


    # check files

    if not os.path.isfile(os.path.join(SUBJECTS_DIR,SUBJECT,"mri","transforms","cc_up.lta")):

        print('WARNING: could not find '+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","transforms","cc_up.lta")+', returning NaNs')

        out = np.empty(N_EIGEN)
        out[:] = np.nan

        return out

    elif not os.path.isfile(os.path.join(SUBJECTS_DIR,SUBJECT,"mri","aseg.mgz")):

        print('WARNING: could not find '+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","aseg.mgz")+', returning NaNs')

        out = np.empty(N_EIGEN)
        out[:] = np.nan

        return out

    elif not os.path.isfile(os.path.join(SUBJECTS_DIR,SUBJECT,"mri","norm.mgz")):

        print('WARNING: could not find '+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","norm.mgz")+', returning NaNs')

        out = np.empty(N_EIGEN)
        out[:] = np.nan

        return out

    ## run make_upright; note: rather than 'make_upright', better use 'mri_cc' 
    ## to compute the transformation matrix, should this ever be necessary.
    #
    #cmd = "make_upright  "+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","norm.mgz")+" "+os.path.join(OUTPUT_DIR,"normCCup.mgz")+" "+os.path.join(OUTPUT_DIR,"cc_up.lta")
    #run_cmd(cmd,"Could not run make_upright")

    # convert lta to xfm (need to adjust some directories when using make_upright)

    cmd = "lta_convert --inlta "+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","transforms","cc_up.lta")+" --outmni "+os.path.join(OUTPUT_DIR,"cc_up.xfm")
    run_cmd(cmd,"Could not convert lta")

    # conduct transform for aseg and norm

    cmd = "mri_convert -i "+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","aseg.mgz")+" -at "+os.path.join(OUTPUT_DIR,"cc_up.xfm")+" -rt nearest -o "+os.path.join(OUTPUT_DIR,"asegCCup.mgz")
    run_cmd(cmd,"Could not conduct cc_up.xfm transform")

    # when using 'make_upright', conducting the transform for norm.mgz is no
    # longer necessary (and will produce the same results)
    cmd = "mri_convert -i "+os.path.join(SUBJECTS_DIR,SUBJECT,"mri","norm.mgz")+" -at "+os.path.join(OUTPUT_DIR,"cc_up.xfm")+" -rt cubic -o "+os.path.join(OUTPUT_DIR,"normCCup.mgz")
    run_cmd(cmd,"Could not conduct cc_up.xfm transform")

    # create fornix mask and surface

    cmd = "mri_binarize --i "+os.path.join(OUTPUT_DIR,"asegCCup.mgz")+" --match 251 252 253 254 255 --dilate 2 --erode 2 --surf "+os.path.join(OUTPUT_DIR,"cc.surf")+" --surf-smooth 1 --o "+os.path.join(OUTPUT_DIR,"cc.mgz")
    run_cmd(cmd,"Could not create fornix mask and surface")

    # --------------------------------------------------------------------------
    # create screenshot

    if CREATE_SCREENSHOT is True:
        createScreenshots(SUBJECT = SUBJECT, SUBJECTS_DIR = SUBJECTS_DIR, 
            INTERACTIVE = False, VIEWS = [('x', -2), ('x', 0), ('x', 2)], LAYOUT = (1, 3),
            BASE = [os.path.join(OUTPUT_DIR,"normCCup.mgz")], OVERLAY = [os.path.join(OUTPUT_DIR,"cc.mgz")], SURF = [os.path.join(OUTPUT_DIR,"cc.surf")], OUTFILE = os.path.join(OUTPUT_DIR,"cc.png"))

    # --------------------------------------------------------------------------
    # run shapeDNA

    if RUN_SHAPEDNA is True:

        import nibabel as nb

        from brainPrintPython import laplaceTria

        surf = nb.freesurfer.io.read_geometry(os.path.join(OUTPUT_DIR,"cc.surf"), read_metadata=True)

        ev, evec = laplaceTria(surf[0], surf[1], k=N_EIGEN)

        d = dict()
        d['Refine'] = 0
        d['Degree'] = 1
        d['Dimension'] = 2
        d['Elements'] = len(surf[1])
        d['DoF'] = len(surf[0])
        d['NumEW'] = N_EIGEN
        d['Eigenvalues'] = ev
        d['Eigenvectors'] = evec

        # return
        return d['Eigenvalues']

    else:

        out = np.empty(N_EIGEN)
        out[:] = np.nan

        return out


