#############################################################################
# Name:    FreeSurferEnv.csh
# Purpose: sets up the environment to run FreeSurfer and FS-FAST 
# Usage:   
#   1. Create an environment variable called FREESURFER_HOME and set it
#      to the directory in which FreeSurfer is installed.
#   2. From a csh or tcsh shell or (.login): 
#      source FreeSurfeEnv.csh
#   3. There are environment variables that should point to locations
#      of software or data used by FreeSurfer. If set prior
#      to sourcing, they will not be changed, but will otherwise be
#      set to default locations:
#        FSFAST_HOME
#        SUBJECTS_DIR
#        MINC_BIN_DIR  
#        MINC_LIB_DIR  
#        FSL_DIR
#   4. If NO_MINC is set (to anything), then all the MINC stuff is ignored.
#   5. If NO_FSFAST is set (to anything), then the startup.m stuff is
#      ignored
#   6. The script will print the final settings for the above variables
#      as well as any warnings about missing directories. If 
#      FS_FREESURFERENV_NO_OUTPUT is set, then no normal output will
#      be made (only error messages).
#
#   The most convenient way to use this script is to write another
#   script that sets FREESURFER_HOME and possibly SUBJECTS_DIR for
#   your set-up, as well as NO_MINC, NO_FSFAST, or
#   FS_FREESURFERENV_NO_OUTPUT as appropriate, and then source this
#   script.
#
#
# $Id: FreeSurferEnv.csh,v 1.9 2005/05/05 18:03:36 kteich Exp $
#############################################################################

set VERSION = '$Id: FreeSurferEnv.csh,v 1.9 2005/05/05 18:03:36 kteich Exp $'

## Get the name of the operating system
set os = `uname -s`
setenv OS $os


## Set this environment variable to suppress the output.
if( $?FS_FREESURFERENV_NO_OUTPUT ) then
    set output = 0
else
    set output = 1
endif
if($?USER == 0 || $?prompt == 0) then
    set output = 0
endif


if( $output ) then
    echo "Setting up enviroment for FreeSurfer/FS-FAST"
    echo $VERSION
endif


## Check if FREESURFER_HOME variable exists, then check if the actual
## directory exists.
if(! $?FREESURFER_HOME) then
  echo "ERROR: environment variable FREESURFER_HOME is not defined"
  echo "       Run the command 'setenv FREESURFER_HOME  FreeSurferHome'"
  echo "       where FreeSurferHome is the directory where FreeSurfer is"
  echo "       installed."
  exit 1;
endif

if(! -e $FREESURFER_HOME) then
  echo "ERROR: $FREESURFER_HOME "
  echo "       does not exist. Check that this value is correct.";
  exit 1;
endif

## Now we'll set directory locations based on FREESURFER_HOME for use
## by other programs and scripts.

## If FS_OVERRIDE is set, this script will automatically assign
## defaults to all locations. Otherwise, it will only do so if the
## variable isn't already set
if(! $?FS_OVERRIDE) then
    setenv FS_OVERRIDE 0
endif


if(! $?FSFAST_HOME || $FS_OVERRIDE) then
  setenv FSFAST_HOME $FREESURFER_HOME/fsfast
endif

if(! $?SUBJECTS_DIR  || $FS_OVERRIDE) then
  setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif

if(! $?NO_MINC && (! $?MINC_BIN_DIR  || $FS_OVERRIDE)) then
  setenv MINC_BIN_DIR $FREESURFER_HOME/minc/bin
endif

if(! $?NO_MINC && (! $?MINC_LIB_DIR  || $FS_OVERRIDE)) then
  setenv MINC_LIB_DIR $FREESURFER_HOME/minc/lib
endif

if(! $?FSL_DIR  || $FS_OVERRIDE) then
  setenv FSL_DIR $FREESURFER_HOME/fsl
endif

setenv FREESURFER_HOME          $FREESURFER_HOME 
setenv FREESURFER_HOME        $FREESURFER_HOME
setenv LOCAL_DIR        $FREESURFER_HOME/local
setenv FUNCTIONALS_DIR  $FREESURFER_HOME/sessions

## Make sure these directories exist.
foreach d ($FSFAST_HOME $SUBJECTS_DIR)
  if(! -e $d ) then
      if( $output ) then
	  echo "WARNING: $d does not exist"
      endif
  endif
end

if( $output ) then
    echo "FREESURFER_HOME $FREESURFER_HOME"
    echo "FSFAST_HOME     $FSFAST_HOME"
    echo "FSL_DIR         $FSL_DIR"
    echo "SUBJECTS_DIR    $SUBJECTS_DIR"
endif

## Talairach subject in anatomical database.
setenv FS_TALAIRACH_SUBJECT talairach


######## --------- Functional Analysis Stuff ----------- #######
if( ! $?NO_FSFAST) then
  setenv FMRI_ANALYSIS_DIR $FSFAST_HOME # backwards compatability
  set SUF = ~/matlab/startup.m
  if(! -e $SUF) then
    echo "INFO: $SUF does not exist ... creating"
    mkdir -p ~/matlab
    touch $SUF

    echo "%------------ FreeSurfer FAST ------------------------%" >> $SUF
    echo "fsfasthome = getenv('FSFAST_HOME');"                     >> $SUF
    echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"       >> $SUF
    echo "path(path,fsfasttoolbox);"                               >> $SUF
    echo "clear fsfasthome fsfasttoolbox;"                         >> $SUF
    echo "%-----------------------------------------------------%" >> $SUF
  endif

  set tmp1 = `grep FSFAST_HOME $SUF       | wc -l`;
  set tmp2 = `grep FMRI_ANALYSIS_DIR $SUF | wc -l`;
  
  if($tmp1 == 0 && $tmp2 == 0) then
      if( $output ) then
	  echo ""
	  echo "WARNING: The $SUF file does not appear to be";
	  echo "         configured correctly. You may not be able"
	  echo "         to run the FS-FAST programs";
	  echo "Try adding the following three lines to $SUF"
	  echo "----------------cut-----------------------"
	  echo "fsfasthome = getenv('FSFAST_HOME');"         
	  echo "fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);"
	  echo "path(path,fsfasttoolbox);"                        
	  echo "clear fsfasthome fsfasttoolbox;"
	  echo "----------------cut-----------------------"
	  echo ""
       endif
  endif
endif

### ----------- MINC Stuff -------------- ####
if(! $?NO_MINC) then
  if(! -d $MINC_BIN_DIR) then
      if( $output ) then
	  echo "WARNING: $MINC_BIN_DIR does not exist.";
      endif
  endif
  if(! -d $MINC_LIB_DIR) then
      if( $output ) then
	  echo "WARNING: $MINC_LIB_DIR does not exist.";
      endif
  endif
  ## Set Load library path ##
  if(! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH  $MINC_LIB_DIR
  else
    setenv LD_LIBRARY_PATH  "$LD_LIBRARY_PATH":"$MINC_LIB_DIR"
  endif
endif

### ----------- FSL ------------ ####
setenv FSL_BIN $FSL_DIR/bin
if(! -d $FSL_BIN) then
    if( $output ) then
	echo "WARNING: $FSL_BIN does not exist.";
    endif
endif


## Set up the path. They should probably already have one, but set a
## basic one just in case they don't. Then add one with all the
## directories we just set.
if(! $?path ) then
  set path = ( ~/bin /bin /usr/bin /usr/local/bin )
endif

set path = ( $FSFAST_HOME/bin     \
             $FREESURFER_HOME/bin/noarch      \
             $FREESURFER_HOME/bin/         \
             $FSL_BIN                   \
	     $path \
            )

if(! $?NO_MINC) then
  set path = ( $MINC_BIN_DIR $path )
endif
rehash;


## Add path to OS-specific dynamic libraries.
if(! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH  $FREESURFER_HOME/lib/
else
    setenv LD_LIBRARY_PATH  "$LD_LIBRARY_PATH":"$FREESURFER_HOME/lib/"
endif


exit 0;
####################################################################
