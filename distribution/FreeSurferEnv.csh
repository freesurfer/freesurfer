#############################################################################
# Name:    FreeSurfeEnv.csh
# Purpose: sets up the environment to run FreeSurfer and FS-FAST 
# Usage:   
#   1. Create an environment variable called FREESURFER_HOME and set it
#      to the directory in which FreeSurfer is installed.
#   2. From a csh or tcsh shell or (.login): 
#      source FreeSurfeEnv.csh
#   3. There are four environment variables that, if set prior
#      to sourcing, will not be changed:
#        FSFAST_HOME
#        FSL_DIR
#        AFNI_DIR  
#        SUBJECTS_DIR
#        NO_MINC
#        MINC_BIN_DIR  
#        MINC_LIB_DIR  
#        NO_FSFAST
#   4. If NO_MINC is set (to anything), then all the MINC stuff is ignored.
#   5. If NO_FSFAST is set (to anything), then the startup.m stuff is
#      ignored
#
# $Id: FreeSurferEnv.csh,v 1.1 2003/05/07 21:52:04 kteich Exp $
#############################################################################

setenv os `uname -s`

## Turn on verboseness if desired ##
if($?FS_VERBOSE) then
  set echo
  set verbose
endif

set VERSION = '$Id: FreeSurferEnv.csh,v 1.1 2003/05/07 21:52:04 kteich Exp $'

echo "Setting up enviroment for FreeSurfer/FS-FAST"
echo $VERSION

## ------- Set up a working path ---------- ##
if(! $?path ) then
  set path = ( . ~/bin /bin /usr/bin /usr/local/bin /usr/bsd \
               /usr/sbin /sbin /usr/etc /usr/ucb)
  rehash;
endif

### ------- Set and check the FREESURFER_HOME directory ----------- ###
if(! $?FREESURFER_HOME) then
  echo "ERROR: environment variable FREESURFER_HOME is not defined"
  echo "Run the command 'setenv FREESURFER_HOME  FreeSurferHome'"
  echo "where FreeSurferHome is the directory where FreeSurfer is installed"
endif
if(! -e $FREESURFER_HOME) then
  echo "ERROR: $FREESURFER_HOME does not exist.  Environment not reset.";
  exit 1;
endif

# ----- Override pre-existing settings -------- #
if(! $?FS_OVERRIDE) setenv FS_OVERRIDE 0;

##------------------------------------------------------------##
if(! $?FSFAST_HOME || $FS_OVERRIDE) then
  setenv FSFAST_HOME $FREESURFER_HOME/fsfast
endif
##------------------------------------------------------------##
if(! $?SUBJECTS_DIR  || $FS_OVERRIDE) then
  setenv SUBJECTS_DIR $FREESURFER_HOME/subjects
endif
##------------------------------------------------------------##
if(! $?NO_MINC && (! $?MINC_BIN_DIR  || $FS_OVERRIDE)) then
  setenv MINC_BIN_DIR $FREESURFER_HOME/minc/bin
endif
##------------------------------------------------------------##
if(! $?NO_MINC && (! $?MINC_LIB_DIR  || $FS_OVERRIDE)) then
  setenv MINC_LIB_DIR $FREESURFER_HOME/minc/lib
endif
##------------------------------------------------------------##
if(! $?AFNI_DIR  || $FS_OVERRIDE) then
  setenv AFNI_DIR $FREESURFER_HOME/afni
endif

##------------------------------------------------------------##
if(! $?FSL_DIR  || $FS_OVERRIDE) then
  setenv FSL_DIR $FREESURFER_HOME/fsl
endif

foreach d ($FSFAST_HOME $SUBJECTS_DIR $AFNI_DIR)
  if(! -e $d )   echo "WARNING: $d does not exist"
end

echo "FREESURFER_HOME $FREESURFER_HOME"
echo "FSFAST_HOME     $FSFAST_HOME"
echo "AFNI_DIR        $AFNI_DIR"
echo "FSL_DIR         $FSL_DIR"
echo "SUBJECTS_DIR    $SUBJECTS_DIR"

setenv CSURF_DIR $FREESURFER_HOME

## Get the name of the operating system ##
set os = `uname -s`

setenv LOCAL_DIR        $CSURF_DIR/local      # utils
setenv FUNCTIONALS_DIR  $CSURF_DIR/sessions   # demo functional ???

## Talairach subject in anatomical database ##
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

### ----------- MINC Stuff -------------- ####
if(! $?NO_MINC) then
  if(! -d $MINC_BIN_DIR) then
    echo "WARNING: $MINC_BIN_DIR does not exist.";
  endif
  if(! -d $MINC_LIB_DIR) then
    echo "WARNING: $MINC_LIB_DIR does not exist.";
  endif
  ## Set Load library path ##
  if(! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH  $MINC_LIB_DIR
  else
    setenv LD_LIBRARY_PATH  "$LD_LIBRARY_PATH":"$MINC_LIB_DIR"
  endif
endif

### ----------- AFNI ------------ ####
setenv AFNI_BIN $AFNI_DIR/$os
if(! -d $AFNI_BIN) then
  echo "WARNING: $AFNI_BIN does not exist.";
endif

### ----------- FSL ------------ ####
setenv FSL_BIN $FSL_DIR/bin
if(! -d $FSL_BIN) then
  echo "WARNING: $FSL_BIN does not exist.";
endif

setenv LOCAL_DIR     $CSURF_DIR/local

### ---------- TCL/TK: Version ------------- ###
## Version 7.4/4.0 is needed for tksurfer stuff
## Version 8 is needed for tony's browser
## We'll set the version to 7.4/4.0 and put tony's
## browser into a wrapper which sets the version to
## 8.0 before running the browser. See browse-sessions.
#setenv TCL_LIBRARY   $LOCAL_DIR/lib/tcl7.4
#setenv TK_LIBRARY    $LOCAL_DIR/lib/tk4.0
setenv TCL_LIBRARY   $LOCAL_DIR/lib/tcl8.3  # tcl startup
setenv TK_LIBRARY    $LOCAL_DIR/lib/tk8.3   # tk startup

setenv TIX_LIBRARY   $LOCAL_DIR/lib/$os/tix4.1

setenv DIAG          0x4040

set path = ( $path \
             $FSFAST_HOME/bin     \
             $FSFAST_HOME/bin/$os \
             $LOCAL_DIR/bin/$os         \
             $CSURF_DIR/bin/noarch      \
             $CSURF_DIR/bin/$os         \
             $FSL_BIN                   \
             $AFNI_BIN                  \
            )

if(! $?NO_MINC) then
  set path = ( $path $MINC_BIN_DIR)
endif
rehash;

# next to be removed!
setenv MRI_DIR $CSURF_DIR  # back compat (test w/o)

if($FS_OVERRIDE) then
  setenv FS_OVERRIDE = 0;
  echo "INFO: reseting FS_OVERRIDE to 0"
endif

if($?FS_VERBOSE) then
  unset echo ;
  unset verbose;
endif

exit 0;
####################################################################
