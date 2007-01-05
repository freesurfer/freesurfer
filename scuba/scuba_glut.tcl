##
## scuba_glut.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:21:46 $
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

SetFrameViewConfiguration 1 c1
set colID [MakeDataCollection Volume]
SetVolumeCollectionFileName $colID /home/kteich/freesurfer/subjects/bert/mri/T1
set layerID [MakeLayer 2DMRI]
SetLayerLabel $layerID "bert"
Set2DMRILayerVolumeCollection $layerID $colID
set viewID [GetViewIDFromFrameColRow 1 0 0]
set level [GetFirstUnusedDrawLevelInView $viewID]
SetLayerInViewAtLevel $viewID $layerID $level
