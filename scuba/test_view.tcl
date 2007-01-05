##
## test_view.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:53 $
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

LoadVolume "/home/kteich/test_data/anatomical/bert.mgz" 1 0
SetLayerLabel 0 "bert"

Set2DMRILayerVolumeCollection 0 0
SetLayerInAllViewsInFrame 0 0

SetViewInPlane 0 y
