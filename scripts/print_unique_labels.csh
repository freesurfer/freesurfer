#! /bin/tcsh -f

# print_unique_labels
#
# Printing unique labels in the input file
#
# Original Author: Lilla Zollei
# Created: 08-09-2011
# MOdifications: 02-07-13: added option to just output a list of lables and not lables + name into a file

alias grep grep
set PrintHelp = 0;
set labelvol = ()
set LOI = -1
set outputfile = ()
set fscolor = $FREESURFER_HOME/FreeSurferColorLUT.txt
set onlylist = 0

set inputargs = ($argv);
set VERSION = 'print_unique_labels.csh @FS_VERSION@';

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

if !($onlylist) then
  echo "Input file: $labelvol" > $outputfile
  echo "Output file: $outputfile"
endif 

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
if ($onlylist) then
  echo $labels
else
  echo "There are $#labels unique lables in this volume. They are: ($labels)"
endif 

if !($onlylist) then
  foreach LOI ($labels)
    # To get FS structure name:
    set linenumber = `sed -n "/^$LOI / =" $fscolor`
    if ($LOI < 10) then
      set structure  = `sed -n "$linenumber s/$LOI   \([0-9a-zA-Z-]*\).*/ \1/p" < $fscolor`
    else if ($LOI < 100) then
      set structure  = `sed -n "$linenumber s/$LOI  \([0-9a-zA-Z-]*\).*/ \1/p" < $fscolor`
    else 
      set structure  = `sed -n "$linenumber s/$LOI \([0-9a-zA-Z-]*\).*/ \1/p" < $fscolor`
    endif
    echo "$structure (label ${LOI})" >> $outputfile
  end
  # more $outputfile
endif

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

    case "--out":
      if ( $#argv < 1) goto arg1err;
      set outputfile = $argv[1]; shift;
      breaksw

    case "--list":
      set onlylist = 1; 
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

  if(($#outputfile == 0) && ($onlylist == 0)) then
    echo "ERROR: must have an output file specified or the list option specified"
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
  echo "USAGE: print_unique_labels"
  echo ""
  echo "Required Arguments:";
  echo "   --vol labelvol   : label volume to be analyzed"
  echo "  and one of the below"
  echo "   --out outputfile : text file where the results are written"
  echo "   --list : only lists the labels"
  echo ""
  echo "Optional Arguments"
  echo ""
  echo "   --version        : print version and exit"
  echo "   --help           : print help and exit"
  echo ""

  if($PrintHelp) \
  cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'

exit 1;

#---- Everything below here is printed out as part of help -----#

BEGINHELP

Prints the list of unique lables (with structure name) in the input volume.

Required Arguments:

--vol labelvol

Label volume to be analyzed

--out outputfile 

Text file where the results are written

or

--list 

Only list the labels (not their anatomical equivalents)

BUGS:
