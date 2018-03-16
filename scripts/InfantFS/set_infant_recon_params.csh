#!/bin/tcsh -ef

#Emma set:
setenv FSSCRIPTSDIR $FREESURFER_HOME/scripts/InfantFS

# training subjects
source $BABYDEVSCRIPTSDIR/atlassubjects.csh

# segmentation postprocessing
# set version = Num4
# set withparcellationlabels = 0 
set parcellationlabels = (9000 9001 9002 9003 9004 9005 9006 9500 9501 9502 9503 9504 9505 9506)

# for the surface extraction step
set leftlabel   = 255
set leftprefix  = left
set rightlabel  = 127
set rightprefix = right

#####
setenv PATH  ${PATH}:${BABYDEVSCRIPTSDIR}:${FSSCRIPTSDIR}
echo $PATH

setenv CLUSTERRUN 0
if (($HOST == launchpad) || ($HOST == tensor)) then
  setenv CLUSTERRUN 1
endif
