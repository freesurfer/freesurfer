#!/bin/bash -p

# Call configuration script:

if [ "`uname -s`" == "Darwin" ]; then
  source $FREESURFER_HOME/SetUpFreeSurfer.sh
fi
