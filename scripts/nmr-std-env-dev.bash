#
# NAME
#
#	nmr-std-env-dev.bash
#
# SYNOPSIS
#
#	source nmr-std-env-dev.bash
#		-- OR --
#	. nmr-std-env-dev.bash
#
# DESCRIPTION
#
#	Sets up the environment to run MGH-NMR standard MRI processing
#	stream: unpacking, functional, structural, and visualisation.
#
#	This file is basically a straightforward port of Doug's (t)csh script
#	nmr-std-env-dev.bash
#
# VERSION
#
#	$Id: nmr-std-env-dev.bash,v 1.2 2004/09/01 21:11:40 rudolph Exp $
#
# AUTHOR
#
#	Rudolph Pienaar - rudolph@nmr.mgh.harvard.edu 
#
# SEE ALSO
#
#	Doug's nmr-std-env-dev.csh
#
# HISTORY
#
# 20 July 2004
# o Initial design and coding.

# Version data
SELF="nmr-std-env-dev.bash"
VERSION='$Id: nmr-std-env-dev.bash,v 1.2 2004/09/01 21:11:40 rudolph Exp $'

## Turn on verboseness if desired ##
if [ -n $NMR_STD_VERBOSE ] ; then
    echo=1
    verbose=1
fi

if [[ $echo == 1 ]] ; then
    echo "$SELF"
    echo "$VERSION"
fi

### -------  Print a warning if the user is inverse ----------- ###
if [ $(whoami) == "inverse" ] ; then
    echo " "
    echo "WARNING: you are currently logged in as user 'inverse'.  This account"
    echo "will eventually be made unavaible to the general NMR community.  Please"
    echo "try logging in as yourself to process your data."
    echo " "
    return 1
fi

# OS data
export OS=$(uname -s)

# Main ENV variables
export FREESURFER_HOME=/space/lyon/1/fsdev/freesurfer_dev
export SUBJECTS_DIR=/cmas/fs/1/users/freesurfer/Subjects
export FSFAST_HOME=/space/lyon/1/fsdev/freesurfer_dev/fsfast

# In the original (t)csh script, the following are OS dependent.
# Here, I assume Linux only
export MINC_BIN_DIR=/space/lyon/9/pubsw/Linux2/packages/mni/current/bin
export MINC_LIB_DIR=/space/lyon/9/pubsw/Linux2/packages/mni/current/lib

# Source the FreeSurfer Environment File
FS_ENV_FILE=$FREESURFER_HOME/FreeSurferEnv.bash
if [ ! -e $FS_ENV_FILE ] ; then
    echo "ERROR: cannot find $FS_ENV_FILE"
    return 1;
fi

source $FS_ENV_FILE

this_shell=$(basename $SHELL)
this_hostname=$(hostname) 
this_host=$(basename $this_hostname .nmr.mgh.harvard.edu)

### Prompt settings. Note that this will only have a meaning if used 
###	in conjunction with Rudolph's ~/.bashrc
if [ $this_shell == bash ] ; then
    export PROMPTPREFIX="[nse]"
fi

if [ -z $MRI_UMASK ] ; then
    export MRI_UMASK=2
fi

return 0;








































