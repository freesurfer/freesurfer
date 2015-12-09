#!/bin/tcsh -ef

# Call configuration script:

if ("`uname -s`" == "Darwin") then
  source $FREESURFER_HOME/SetUpFreeSurfer.csh
endif
