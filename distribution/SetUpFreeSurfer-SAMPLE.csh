# SetUpFreeSurfer.csh
#

# This is a sample SetUpFreeSurfer.csh file. You should be able to use
# it with just a few modifications.

# Set this to the location of freesrfer/.
if (! $?FREESURFER_HOME) then 
    setenv FREESURFER_HOME REPLACE_WITH_PREFIX
endif    
 
# Set this to your subjects/ dir, usually freesurfer/subjects/
if (! $?SUBJECTS_DIR) then
    setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif

# Set this to your functional sessions dir, usually freesurfer/sessions/
if (! $?FUNCTIONALS_DIR) then
    setenv FUNCTIONALS_DIR $FREESURFER_HOME/sessions
endif

# Specify the location of the MINC tools, such as...
#setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
#setenv MINC_LIB_DIR /usr/pubsw/packages/mni/current/lib
# ... or disable them.
setenv NO_MINC

# Enable or disable fsfast.
setenv NO_FSFAST

# Call configuration script.
source $FREESURFER_HOME/FreeSurferEnv.csh

