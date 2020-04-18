#!/bin/tcsh -ef

# training subjects
source $FREESURFER_HOME/bin/atlassubjects.csh
# parcellation labels
set parcellationlabels = (9000 9001 9002 9003 9004 9005 9006 9500 9501 9502 9503 9504 9505 9506)

# for the surface extraction step
set leftlabel   = 255
set leftprefix  = left
set rightlabel  = 127
set rightprefix = right

set regtype     = NIFTYREG;
#####

# Note, the below should be personalized if run outside of the Martinos
setenv CLUSTERRUN 0
if (($HOST == launchpad) || ($HOST == tensor)) then
  setenv CLUSTERRUN 1
endif
