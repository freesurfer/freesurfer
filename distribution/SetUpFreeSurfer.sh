#!/bin/bash -p

#
# SetUpFreeSurfer.sh
#

# This is a sample SetUpFreeSurfer.sh file.
# Edit as needed for your specific setup.
# The defaults should work with most installations.

# Set this to the location of the freesurfer installation.
if [ -z $FREESURFER_HOME ]; then
    echo " ERROR: Environment variable FREESURFER_HOME must be defined prior to sourcing Freesurfer." 
    return
fi

# Set this to your subjects/ dir, usually freesurfer/subjects/
if [ -z $SUBJECTS_DIR ]; then
    export SUBJECTS_DIR=$FREESURFER_HOME/subjects
fi

# Set this to your functional sessions dir, usually freesurfer/sessions/
if [ -z $FUNCTIONALS_DIR ]; then
    export FUNCTIONALS_DIR=$FREESURFER_HOME/sessions
fi

# Specify the location of the MINC tools.
# Necessary only if the script $FREESURFER_HOME/FreeSurferEnv.csh
# does not find the tools (and issues warnings pertaining to
# the following two environment variables, which have example
# locations that might need user-specific modification):
#export MINC_BIN_DIR=/usr/pubsw/packages/mni/current/bin
#export MINC_LIB_DIR=/usr/pubsw/packages/mni/current/lib
# ... or just disable the MINC toolkit (although some Freesurfer
# utilities will fail!)
#export NO_MINC=1

# Enable or disable fsfast (enabled by default)
#export NO_FSFAST=1

# set user prompt to optional input string
if [ $# -eq 1 ]; then
  PS1="(fsenv-$1) [\h:\W]"
fi

# Call configuration script:
source $FREESURFER_HOME/FreeSurferEnv.sh
