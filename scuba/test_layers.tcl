##
## test_layers.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:51 $
##    $Revision: 1.2 $
##
## Copyright (C) 2002-2007,
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
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
