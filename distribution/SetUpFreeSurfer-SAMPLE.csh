# SetUpFreeSurfer.csh
#

# This is a sample SetUpFreeSurfer.csh file. You should be able to use
# it with just a few modiciations.

# Set this to the location of freesrfer/.
if (! $?FREESRFER_HOME) then 
    setenv FREESURFER_HOME REPLACE_WITH_PREFIX
endif    
 
# Call configuration script.
pushd $FREESURFER_HOME
source FreeSurferEnv.csh
popd

# Set this to your subjects/ dir, usually freesurfer/subjects/
if (! $?SUBJECTS_DIR) then
    setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif

# Specify the location of the MINC tools...
#setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
#setenv MINC_LIB_DIR /usr/pubsw/packages/mni/current/lib
# ... or disable them.
setenv NO_MINC

# Enable or disable fsfast.
setenv NO_FSFAST
