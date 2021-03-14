##
## twocond-masked-views.tcl
##
##
## Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

if [info exists noexit] { unsetenv condition0 }
set mask ${hemi}-cortex
source $env(FREESURFER_HOME)/lib/tcl/twocond-views.tcl
#mask_label ${hemi}-cortex
