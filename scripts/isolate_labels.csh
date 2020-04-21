#! /bin/tcsh -f

# isolate_labels
#
# Writing out all / some labels in a binary file format
# (for subsequent shape analysis)
#
# Original Author: Lilla Zollei
# Created: 02-14-2011


alias grep grep
set doAll = 1;
set PrintHelp = 0;
set labelvol = ()
set LOI = -1
set outprefix = ()
set keepvalues = 0

set inputargs = ($argv);
set VERSION = 'isolate_labels.csh @FS_VERSION@';

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

echo "Input file: $labelvol"
echo "Output prefix: $outprefix"

# Working with some fsl code so make sure that there is a NIFTI
# version of the input
set didconversion = 0
if ( ${labelvol:e} != gz) then
 set cmd = (mri_convert $labelvol ${labelvol:r}.nii.gz)
 # echo $cmd
 $cmd
 set labelvol = ${labelvol:r}.nii.gz
 set didconversion = 1
endif

if ($doAll) then # Find all the labels in the volume
  set labels = ();
  set tmp = `fslstats $labelvol -R`
  set min  = ${tmp[1]:r}
  set max  = ${tmp[2]:r} 
  @ bins = $max - $min # each label gets a bin
  @ bins = $bins + 1
  set histo = `fslstats $labelvol -H $bins $min $max`
  # echo $histo
  # Now find all the non-zero elements in the histogram
  set counter = 1
  while ($counter < $bins + 1) 
    if (${histo[$counter]:r} > 0) then
      @ l = $counter - 1
      set labels = ($labels $l)
    endif
    @ counter = $counter + 1
  end
  echo "There are $#labels unique lables in this volume. They are: ($labels)"
else
  echo "Label of interest: $LOI"
  set labels = $LOI
endif

##
## mri_binarize --i --match --o 
##

foreach LOI ($labels)

  @ lower = $LOI - 1
  @ upper = $LOI + 1

  set tmp = `fslstats $labelvol -l $lower -u $upper -v`
  set labelN   = ${tmp[1]:r}
  set labelVol = $tmp[2]
  echo "Number of voxels and volume for label ${LOI} is = (${labelN}, ${labelVol}mm^3)"

  if ($keepvalues) then
    set cmd = (mri_binarize --i $labelvol --match $LOI --o ${outprefix}_label$LOI.mgz --binval $LOI)
  else
    set cmd = (mri_binarize --i $labelvol --match $LOI --o ${outprefix}_label$LOI.mgz )
  endif
  echo $cmd
  eval $cmd

end

if ($didconversion) then
  rm $labelvol
endif

exit 0;

############--------------##################
############--------------##################
parse_args:
set cmdline = ($argv);
while( $#argv != 0 )

  set flag = $argv[1]; shift;

  switch($flag)

    case "--vol":
      if ( $#argv < 1) goto arg1err;
      set labelvol = $argv[1]; shift;
      breaksw

    case "--outprefix":
      if ( $#argv < 1) goto arg1err;
      set outprefix = $argv[1]; shift;
      breaksw

    case "--keepval":
      set keepvalues = 1; 
      breaksw

    case "--L":
    case "--l":
      if ( $#argv < 1) goto arg1err;
      set LOI = $argv[1]; shift;
      set doAll = 0;
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

  if($#labelvol == 0) then
    echo "ERROR: must have a label volume to analyze"
    exit 1;
  endif

  if($#outprefix == 0) then
    echo "ERROR: must have an output file prefix specified"
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
  echo "USAGE: isolate_labels"
  echo ""
  echo "Required Arguments:";
  echo "   --vol labelvol        : label volume to be analyzed"
  echo "   --outprefix outprefix : prefix of binary label file(s)"
  echo ""
  echo "Optional Arguments"
  echo ""
  echo "   --L label        : the particular label to be analyzed"
  echo "   --l label        : the particular label to be analyzed"
  echo "   --version        : print version and exit"
  echo "   --keepval        : keeps original label values"
  echo "   --help           : print help and exit"
  echo ""

  if($PrintHelp) \
  cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1;

#---- Everything below here is printed out as part of help -----#

BEGINHELP

Separates out a particular or all lables into individual binary files 
(for subsequent shape analysis).


Required Arguments:

--vol labelvol

Label volume to be analyzed

--outprefix outprefix

Text file where the results are written

Optional Arguments:

--L label
--l label

The particular label to be analyzed. By default, it is ALL labels in the volume.

BUGS:
