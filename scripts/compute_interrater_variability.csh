#! /bin/tcsh -f

# compute_interrater_variability
#
# Computing interrater variability.
#
# Original Author: Lilla Zollei
# Created: 04-27-2010

alias grep grep
set doAll = 1;
set PrintHelp = 0;
set labelvol = ()
set LOI = -1
set outputprefix = ()

set inputargs = ($argv);
set VERSION = 'compute_interrater_variability.csh @FS_VERSION@';

if($#argv == 0) goto usage_exit;
set n = `echo $argv | egrep -e --version | wc -l`
if($n != 0) then
  echo $VERSION
  exit 0;
endif
set n = `echo $argv | egrep -e --help | wc -l`
if($n != 0) then
  set PrintHelp = 1;
  goto usage_exit;
endif

source $FREESURFER_HOME/sources.csh

goto parse_args;
parse_args_return:
goto check_params;
check_params_return:

############--------------##################
set DateString = "`date '+%y%m%d%H%M'`"

set LF = mri_interrater_variability.$DateString.log
if(-e $LF) mv $LF $LF.old;
echo ""
echo "Log file is $LF"
echo ""

echo "Logfile for compute_interrater_variability" >> $LF
date |& tee -a $LF
echo $inputargs |& tee -a $LF
echo $VERSION |& tee -a $LF
hostname |& tee -a $LF
uname -a |& tee -a $LF

############--------------##################

set meanhaus = ${outputprefix}.mean.haus.txt
set cmd = (mri_hausdorff_dist $file1 $file2 $meanhaus)
echo $cmd |& tee -a $LF
$cmd |& tee -a $LF
if($status) exit 1;

set maxhaus = ${outputprefix}.max.haus.txt
set cmd = (mri_hausdorff_dist -max $file1 $file2 $maxhaus)
echo $cmd |& tee -a $LF
$cmd |& tee -a $LF
if($status) exit 1;

set overlapoutput = ${outputprefix}.overlap.txt
set cmd = (mri_compute_overlap -a -s -l $overlapoutput $file1 $file2)
echo $cmd |& tee -a $LF
$cmd |& tee -a $LF
if($status) exit 1;

exit 0;

############--------------##################
############--------------##################
parse_args:
set cmdline = ($argv);
while( $#argv != 0 )

  set flag = $argv[1]; shift;

  switch($flag)

    case "--vol1":
      if ( $#argv < 1) goto arg1err;
      set file1 = $argv[1]; shift;
      breaksw

    case "--vol2":
      if ( $#argv < 1) goto arg1err;
      set file2 = $argv[1]; shift;
      breaksw

    case "--out":
      if ( $#argv < 1) goto arg1err;
      set outputprefix = $argv[1]; shift;
      breaksw

    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw

end

goto parse_args_return;
############--------------##################
############--------------##################
check_params:

  if(($#file1 == 0) || ($#file2 == 0)) then
    echo "ERROR: must have two label volumes indicated as an input"
    exit 1;
  endif

  if($#outputprefix == 0) then
    echo "ERROR: must have an output prefix file named"
    exit 1;
  endif

goto check_params_return;
############--------------##################
############--------------##################
arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1
############--------------##################
############--------------##################
usage_exit:
  echo "USAGE: compute_interrater_variability"
  echo ""
  echo "Required Arguments:";
  echo "   --vol1 labelvol1   : rater1 label volume"
  echo "   --vol1 labelvol2   : rater2 label volume"
  echo "   --out outputprefix : text file prefix where the results are written; a total of 3 files are produced"
  echo "   --version          : print version and exit"
  echo "   --help             : print help and exit"
  echo ""

  if($PrintHelp) \
  cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1;

#---- Everything below here is printed out as part of help -----#

BEGINHELP

Computes the interrater variability between the input files. The input volumes are 
expected to be label files segmented by two different users or the same user at 
different time points. The three types of comparisons that are made are the following:
(1) mean hausdorff distance (smaller number means higher consistency)
(2) max hausdorff distance  (smaller number means higher consistency)
(3) label volume difference, Dice and Jaccard overlap measures (smaller difference and higher overlap 
    means higher consistency)

Required Arguments:

--vol1 labelvol1
--vol2 labelvol2

Label volumes to be analyzed

--out outputprefix 

Text file prefix where the results are written. A total of three files are produced.

BUGS:
