#! /bin/tcsh -f

# mri_reorient_LR
#
# Script to reorient (baby) brain scans for easier segmentation.
#
# Original Author: Lilla Zollei
# Created: 07-16-2010

set inputargs = ($argv);
set VERSION = 'mri_reorient_LR.csh @FS_VERSION@';

set inputvol      = ();
set outputvol     = ();
set outreg        = ();
set displayresult = 1;
set newinputnii   = 0;
set newoutputnii  = 0;
set cleanup       = 0;
set PrintHelp     = 0;

# Parsing and checking the input arguments
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
############--------------##################

#set orientation = `mri_info --orientation $inputvol`
#if ($orientation == RIA) then
#  set neworientation = LIA
#else
#  set neworientation = RIA
#endif

#switch ($orientation)
# case 'RIA':
#   set neworientation = LIA;
#   breaksw;
# case 'LIA':
#   set neworientation = RIA;
#   breaksw;
# case 'RAS':
#   set neworientation = LAS;
#   breaksw;
# case 'LAS':
#   set neworientation = RAS;
#   breaksw;
# default:
#   echo "***Input orientation ($orientation) not handled! Sorry..."
#   exit 1;
#endsw

#echo "registration will take place between: $orientation and $neworientation"

# If input does not have nii(.gz) then do conversion -- FLIRT canot handle .mgz
echo "***Input check..."
set originputvol = $inputvol
if ( (${inputvol:e} != "nii") && (! (${inputvol:e} == "gz" && ${inputvol:r:e} == "nii") ) ) then
   set cmd = ( mri_convert $inputvol ${inputvol:r}.nii.gz )
   $cmd
   set inputvol    = ${inputvol:r}.nii.gz
   set newinputnii = 1
endif

### FLIP input volume around lh/rh axis
echo "***Flip..."
set outputdir = ${outputvol:h}
#set inputvolLR = $outputdir/${inputvol:t:r:r}.lrflipped.nii.gz
set inputvolLR = $outputdir/${inputvol:t:r:r}.lrreversed.nii.gz
#set cmd = (mri_convert --out_orientation $neworientation $originputvol $inputvolLR)
set cmd = (mri_convert --left-right-reverse $originputvol $inputvolLR)
echo $cmd; eval $cmd

### REGISTER between orig and lhrhflipped volumes
echo "***Register..."
set flirtoutput = ${inputvolLR:r:r}.flirt6DOFreg.nii.gz
set flirtfslmat = ${flirtoutput:r:r}.fslmat
set cmd = (flirt.fsl -dof 6 -in $inputvol -ref $inputvolLR -out $flirtoutput -omat $flirtfslmat)
echo $cmd; eval $cmd

### COMPUTE half angle rotation
echo "***Half angle..."
set avscale = `avscale $flirtfslmat`
echo $avscale > $outputdir/avscale.txt
set matrix = `perl -0777 -ne 'print $1 . "\n" while s/Forward half transform = (.*?)Backward//s' $outputdir/avscale.txt`
set newtransformation = ${flirtoutput:r:r}.HALFFORW.fslmat 
echo $matrix[1]  $matrix[2]  $matrix[3]  $matrix[4]  >  $newtransformation
echo $matrix[5]  $matrix[6]  $matrix[7]  $matrix[8]  >> $newtransformation
echo $matrix[9]  $matrix[10] $matrix[11] $matrix[12] >> $newtransformation
echo $matrix[13] $matrix[14] $matrix[15] $matrix[16] >> $newtransformation

### APPLYING half-angle rotation
echo "***Create output..."
if (( ${outputvol:e} != "nii") && (! (${outputvol:e} == "gz" && ${outputvol:r:e} == "nii") )) then
  set origoutputvol = $outputvol
  set outputvol = ${outputvol:r}.nii.gz
  set newoutputnii = 1
endif
set cmd = (flirt -applyxfm -in $inputvol -ref $inputvolLR -out $outputvol -init $newtransformation)
echo $cmd; eval $cmd
if ($newoutputnii) then
  set cmd = (mri_convert $outputvol $origoutputvol)
  $cmd
  set outputvol = $origoutputvol
endif

### SAVE
if ($#outreg > 0) then
  if (${outreg:e} == lta) then
    set cmd = (tkregister2_cmdl --mov $inputvol --targ $inputvolLR --ltaout $outreg --fsl $newtransformation --reg bogus.reg)
  else # assume fslmat
    set cmd = (cp $newtransformation $outreg)
  endif
  echo $cmd; eval $cmd
endif

### DISPLAY
if ($displayresult) then 
  freeview -v  $originputvol $outputvol
endif

### CLEANING UP
if ($cleanup) then 
  set cmd = (rm -f $inputvolLR $flirtoutput $flirtfslmat $outputdir/avscale.txt $newtransformation )
  echo $cmd; eval $cmd
endif 
if ($newinputnii) then
  set cmd = (rm -f $inputvol)
  echo $cmd; eval $cmd
endif
if ($newoutputnii) then
  set cmd = (rm -f ${outputvol:r}.nii.gz)
  echo $cmd; eval $cmd
endif

exit 0;

##############

############--------------##################
############--------------##################
parse_args:
set cmdline = ($argv);
while( $#argv != 0 )

  set flag = $argv[1]; shift;

  switch($flag)

    case "--i":
      if ( $#argv < 1) goto arg1err;
      set inputvol = $argv[1]; shift;
      breaksw

    case "--o":
      if ( $#argv < 1) goto arg1err;
      set outputvol = $argv[1]; shift;
      breaksw

    case "--outreg":
      if ( $#argv < 1) goto arg1err;
      set outreg = $argv[1]; shift;
      breaksw

    case "--disp":
      if ( $#argv < 1) goto arg1err;
      set displayresult = $argv[1]; shift;
      breaksw

    case "--clean":
      if ( $#argv < 1) goto arg1err;
      set cleanup = $argv[1]; shift;
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

  if($#inputvol == 0) then
    echo "ERROR: must indicate an input vol"
    exit 1;
  endif

  if($#outputvol == 0) then
    echo "ERROR: must indcate an output vol"
    exit 1;
  endif

goto check_params_return;
############--------------##################


############--------------##################
usage_exit:
  echo "USAGE: mri_reorient_LR"
  echo ""
  echo "Required Arguments:";
  echo "   --i vol : input file to be reoriented"
  echo "   --o vol : reoriented input file"
  echo "Optional Arguments"
  echo ""
  echo "   --disp       : display registration result using FreeView (def = 1)"
  echo "   --clean      : delete all aux and reg files (def = 0)"
  echo "   --outreg     : write out the registration file that is applied to the reoriented input file (fslmat or lta)"
  echo "   --version    : print version and exit"
  echo "   --help       : print help and exit"
  echo ""

  if($PrintHelp) \
  cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1;
