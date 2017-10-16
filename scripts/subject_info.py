# Original author - Krish Subramaniam
# $Id: subject_info.py,v 1.3 2009/04/11 19:39:14 krish Exp $

from __future__ import print_function
import os, sys

__all__ = ['check_subjdirs']

def check_subjdirs():
    """
    Quit if SUBJECTS_DIR is not defined as an environment variable. This is not
    a function which returns a boolean. Execution is stopped if not found.
    If found, returns the SUBJECTS_DIR
    """
    if 'SUBJECTS_DIR' not in os.environ:
        print('ERROR: SUBJECTS_DIR environment variable not defined!')
        sys.exit(1)
    return os.environ['SUBJECTS_DIR']
        
