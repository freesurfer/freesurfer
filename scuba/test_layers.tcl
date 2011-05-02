##
## test_layers.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/05/02 21:02:44 $
##    $Revision: 1.3 $
##
## Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

LoadVolume "/home/kteich/test_data/anatomical/25cubed.mgz" 1 0
SetLayerLabel 0 "Layer 0"

for  { set n 1 } { $n < 10 } { incr n } {
    Make2DMRILayer "Layer $n"
    Set2DMRILayerVolumeCollection $n 0
    SetLayerInAllViewsInFrame 0 $n
}
