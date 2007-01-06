#! /bin/tcsh -f

#
# inu_correct_vol.csh
#
# REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
#
# Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR
# CVS Revision Info:
#    $Author: nicks $
#    $Date: 2007/01/06 00:01:14 $
#    $Revision: 1.2 $
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


set vol_in = $1
set vol_out = $2

mri_convert ${vol_in} input.mnc
nu_correct -clobber input.mnc nu_1.mnc
nu_correct -clobber nu_1.mnc nu_2.mnc

mri_convert nu_2.mnc ${vol_out}

rm -f input*
rm -f nu_1.mnc
rm -f *.imp
rm -f nu_2.*
