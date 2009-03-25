# Original author - Krish Subramaniam
# $Id: subject_info.py,v 1.2.2.2 2009/03/25 19:54:33 krish Exp $

import os, sys

__all__ = ['check_subjdirs']

def check_subjdirs():
    """
    Quit if SUBJECTS_DIR is not defined as an environment variable. This is not
    a function which returns a boolean. Execution is stopped if not found.
    If found, returns the SUBJECTS_DIR
    """
    if 'SUBJECTS_DIR' not in os.environ:
        print 'ERROR: SUBJECTS_DIR environment variable not defined!'
        sys.exit()
    return os.environ['SUBJECTS_DIR']
        
