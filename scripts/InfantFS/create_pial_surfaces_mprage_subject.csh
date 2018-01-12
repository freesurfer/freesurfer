#!/bin/tcsh -ef

## read input
set cmdline = ($argv);
while( $#argv != 0 )
  set flag = $argv[1]; shift;
  switch($flag)
    case "--s":
      if ( $#argv < 1) goto arg1err;
      set subj = $argv[1]; shift;
      breaksw
    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end
## check input
if($#subj == 0) then
  echo "ERROR: must spec a subject id"
  exit 1;
endif

######## PROCESSING

echo Processing the following subject: $subj

set workdir = $SUBJECTS_DIR/$subj/work 
set surfdir = $SUBJECTS_DIR/$subj/surf
set vol     = $SUBJECTS_DIR/$subj/mri/aseg.nii.gz
set segfile = $vol 


foreach hemi (lh rh)
  set cmd = ( mris_make_surfaces -grad_dir 1 -intensity .3 -output .tmp -pial_offset .25 -nowhite -noaparc -cover_seg $vol -orig_pial white $subj $hemi)
  echo $cmd; eval $cmd
  ### BRUCE's new command
  #set cmd = ( mris_make_surfaces -pa 32 2 -grad_dir 1 -intensity .3 -output .brf -max_thickness 10 -pial_offset .25 -nowhite -noaparc -cover_seg $vol -orig_pial white $subj $hemi)
  echo $cmd; eval $cmd
  ##
  set cmd = ( mris_smooth -nw -n 2 $surfdir/$hemi.pial.tmp $surfdir/$hemi.sm.pial.tmp )
  echo $cmd; eval $cmd
  set cmd = (mv $surfdir/$hemi.sm.pial.tmp $surfdir/$hemi.pial)
  echo $cmd; eval $cmd
  set cmd = (mv $surfdir/$hemi.thickness.tmp $surfdir/$hemi.thickness)
  echo $cmd; eval $cmd
  ### NOTE: cleanup!!
end
