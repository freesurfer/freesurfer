##
## dara_data.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:03 $
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

LoadFunctionalOverlay /space/beast/1/users/dara/normalsnew/bold/edpanalmc/tal-ffx/allvfix sig /space/beast/1/users/dara/normalsnew/bold/edpanalmc/tal-ffx/register.data
Overlay_SetThreshold 2 5 1
SetCursor 0 96 67 128
SelectValuesByFuncVoxel 1
RedrawScreen
