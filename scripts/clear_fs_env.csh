#! /bin/tcsh -ef

#
# clear_fs_env.csh
#
# Original Author: Nick Schmansky
# CVS Revision Info:
#    $Author: nicks $
#    $Date: 2007/01/06 00:01:13 $
#    $Revision: 1.4 $
#
# Copyright (C) 2002-2007,
# The General Hospital Corporation (Boston, MA).
# All rights reserved.
#
# Distribution, usage and copying of this software is covered under the
# terms found in the License Agreement file named 'COPYING' found in the
# FreeSurfer source code root directory, and duplicated here:
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
#
# General inquiries: freesurfer@nmr.mgh.harvard.edu
# Bug reports: analysis-bugs@nmr.mgh.harvard.edu
#

unsetenv FREESURFER_HOME
unsetenv FSFAST_HOME
unsetenv SUBJECTS_DIR
unsetenv FUNCTIONALS_DIR
unsetenv MINC_BIN_DIR
unsetenv MINC_LIB_DIR
unsetenv MNI_DATAPATH
unsetenv MNI_PERL5LIB
unsetenv PERL5LIB
unsetenv GSL_DIR
unsetenv VXL_DIR
unsetenv QTDIR
unsetenv FS_TCL_LIB_DIR
unsetenv TCLLIBPATH
unsetenv TCL_LIBRARY
unsetenv TK_LIBRARY
unsetenv TIX_LIBRARY
unsetenv BLT_LIBRARY
unsetenv MISC_LIB
unsetenv FSL_DIR
unsetenv FSLDIR
unsetenv LD_LIBRARY_PATH
unsetenv DYLD_LIBRARY_PATH
