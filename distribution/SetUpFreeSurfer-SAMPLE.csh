# SetUpFreeSurfer.csh
#

# This is a sample SetUpFreeSurfer.csh file. You should be able to use
# it with just a few modiciations.

# Set this to the location of freesrfer/.
setenv FREESURFER_HOME <full path to freesurfer/>
     
# Call configuration script.
cd $FREESURFER_HOME
source FreeSurferEnv.csh

# Set this to your subjects/ dir, usually freesurfer/subjects/
setenv SUBJECTS_DIR $FREESURFER_HOME/subjects

# Specify the location of the MINC tools...
#setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
#setenv MINC_LIB_DIR /usr/pubsw/packages/mni/current/lib
# ... or disable them.
setenv NO_MINC

# Enable or disable fsfast.
setenv NO_FSFAST
