#!/bin/tcsh -ef

##
## For infants between 0-2yrs of age!
##

set inputargs = ($argv);
set fsdirset = 0;

## read input
set cmdline = ($argv);
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--s":
      if ( $#argv < 1) goto arg1err;
      set subj = $argv[1]; shift;
      breaksw
    case "--outdir":
      if ( $#argv < 1) goto arg1err;
      set FS_DIR = $argv[1]; shift;
      set fsdirset = 1;
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end

##
## check input
##
if($#SUBJECTS_DIR == 0) then
  echo "ERROR: must spec a SUBJECTS_DIR to indicate the location of the input data"
  exit 1;
endif
if($#subj == 0) then
  echo "ERROR: must spec a subject id"
  exit 1;
endif

echo Processing the following subject: $subj

if !($fsdirset) then
  set FS_DIR = $SUBJECTS_DIR/
endif

set missing = ()

set mridir = $FS_DIR/$subj/mri
foreach f (norm aseg aparc+aseg)
  set file = $mridir/$f.mgz
  if (! -e $file) then
    set missing = ($missing $file)
  endif
end

set surfdir = $FS_DIR/$subj/surf
foreach hemi (lh rh)
  foreach f (inflated pial sphere white smoothwm sulc)
    set file = $surfdir/$hemi.$f
    if (! -e $file) then
      set missing = ($missing $file)
    endif
  end
end

set surfdir = $SUBJECTS_DIR/$subj/surf
foreach hemi (lh rh)
  foreach f (inflated.H inflated.K)
    set file = $surfdir/$hemi.$f
    if ((! -e $file) && (! -e ${file:r}))  then 
      set missing = ($missing $file)
    endif
  end
end

if ($#missing > 0) then
  echo "The following key files are missing: "
  echo $missing
  echo "You need to find these files or run infant_reconall again on the data"
else
  echo "It seems that no key files are missing."
endif

exit 0;
