#! /bin/tcsh -f

#
# inu_correct_vol.csh
#
#
#
# Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
#
# Terms and conditions for use, reproduction, distribution and contribution
# are found in the 'FreeSurfer Software License Agreement' contained
# in the file 'LICENSE' found in the FreeSurfer distribution, and here:
#
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
#
# Reporting: freesurfer@nmr.mgh.harvard.edu
#
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
