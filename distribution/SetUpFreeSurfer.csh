#!/bin/tcsh -ef

#
# SetUpFreeSurfer.csh
#

# This is a sample SetUpFreeSurfer.csh file. 
# Edit as needed for your specific setup.  
# The defaults should work with most installations.

# Set this to the location of the freesurfer installation.
if (! $?FREESURFER_HOME) then 
    echo " ERROR: Environment variable FREESURFER_HOME must be defined prior to sourcing Freesurfer." 
    exit
endif    

if ($?FREESURFER_FSPYTHON) then 
   if ( -d $FREESURFER_FSPYTHON ) then
      set path = ( $FREESURFER_FSPYTHON/bin $path )
   endif
else if ( -d $FREESURFER_HOME/../fspython/bin ) then
   setenv FREESURFER_FSPYTHON $FREESURFER_HOME/../fspython
   set path = ( $FREESURFER_FSPYTHON/bin $path )
endif

# Set this to your subjects/ dir, usually freesurfer/subjects/
if (! $?SUBJECTS_DIR) then
    setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif

# Set this to your functional sessions dir, usually freesurfer/sessions/
if (! $?FUNCTIONALS_DIR) then
    setenv FUNCTIONALS_DIR $FREESURFER_HOME/sessions
endif

# Specify the location of the MINC tools.
# Necessary only if the script $FREESURFER_HOME/FreeSurferEnv.csh
# does not find the tools (and issues warnings pertaining to
# the following two environment variables, which have example
# locations that might need user-specific modification):
#setenv MINC_BIN_DIR /usr/pubsw/packages/mni/current/bin
#setenv MINC_LIB_DIR /usr/pubsw/packages/mni/current/lib
# ... or just disable the MINC toolkit (although some Freesurfer
# utilities will fail!) 
#setenv NO_MINC

# Enable or disable fsfast (enabled by default)
#setenv NO_FSFAST

# Call configuration script:
source $FREESURFER_HOME/FreeSurferEnv.csh
