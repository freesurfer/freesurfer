# fast-log-preamble
# $Id: fast-log-preamble.csh,v 1.2 2015/12/10 22:51:29 zkaufman Exp $
#
# This file is meant to be sourced at the start of one of the FAST
# scripts.  It prints out generic information that should be stored 
# in the a log file.
source $FREESURFER_HOME/sources.csh

if(! $?FMRI_ANALYSIS_DIR) then
  echo "WARNING: FMRI_ANALYSIS_DIR environment variable is not defined"
  set FASTVER = "unknown";
else
  echo "FMRI_ANALYSIS_DIR: $FMRI_ANALYSIS_DIR"
  if(! -e $FMRI_ANALYSIS_DIR) then
    echo "WARNING: FMRI_ANALYSIS_DIR $FMRI_ANALYSIS_DIR does not exist"
    set FASTVER = "unknown";
  else
    set FASTVERFILE = "$FMRI_ANALYSIS_DIR/docs/version";
    if(! -e $FASTVERFILE) then
      echo "WARNING: $FASTVERFILE does not exist"
      set FASTVER = "unknown";
    else
      set FASTVER = `cat $FMRI_ANALYSIS_DIR/docs/version`;
    endif
  endif
endif

echo "FAST Version $FASTVER"  
echo "User: `id`"            
echo "CurrentDir: `pwd`"
echo "Command: $0"     
echo "Arguments: $argv"  
echo "Machine: `uname -a`"
echo "Date: `date`"

if($?SUBJECTS_DIR) then
  echo "SUBJECTS_DIR: $SUBJECTS_DIR"
else
  echo "SUBJECTS_DIR: not defined"
endif

echo "  Comments or questions: analysis-bugs@nmr.mgh.harvard.edu"

#### done ######
