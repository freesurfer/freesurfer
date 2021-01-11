#!/bin/tcsh -ef

#source subjects.csh

set leftlabel   = 255
set leftprefix  = left
set rightlabel  = 127
set rightprefix = right

# processing stages 
set getlabels   = 1 # NEEDS TO BE ON ALL THE TIME! (--> alllables)
set getwmhemis  = 1
set getfilledwm = 1
set surfaceprep = 1
set getwhitesurfaces = 1
#

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
mkdir -p $workdir
set surfdir = $SUBJECTS_DIR/$subj/surf
set mridir = $SUBJECTS_DIR/$subj/mri
set labeldir = $SUBJECTS_DIR/$subj/label 
set vol     = $SUBJECTS_DIR/$subj/mri/aseg.nii.gz
set segfile = $vol 

###
### Get list of all existing labels
###
if ($getlabels) then
  set niisegfile = $segfile 
  set alllabels = ();
  set tmp = `fslstats $niisegfile -R`
  set min  = ${tmp[1]:r}           
  set max  = ${tmp[2]:r}           
  @ bins = $max - $min # each label gets a bin
  @ bins = $bins + 1                          
  set histo = `fslstats $niisegfile -H $bins $min $max`
  # Now find all the non-zero elements in the histogram
  set counter = 1                                      
  while ($counter < $bins + 1)                         
    if (${histo[$counter]:r} > 0) then                 
      @ l = $counter - 1                               
      set alllabels = ($alllabels $l)                        
    endif                                              
    @ counter = $counter + 1                           
  end 
  echo $alllabels
endif #getlabels

###
### construct the WM hemi volumes (127, 255)
###
if ($getwmhemis) then
  # Left
  set filelist = ()
  foreach label ( 2 9 10 11 12 13 25 26 28 30) 
    set cmd = (mri_binarize --i $vol --o $workdir/${vol:t:r:r}.$leftprefix.$label.mgz --min $label --max $label --binval $leftlabel)
    foreach temp ($alllabels) 
      if ($temp == $label) then
        set filelist = ($filelist $workdir/${vol:t:r:r}.$leftprefix.$label.mgz)
        if !(-e $workdir/${vol:t:r:r}.$leftprefix.$label.mgz) then
          echo $cmd; eval $cmd
        endif
        break
      endif
    end
  end

  set cmd = (mri_or $filelist $workdir/${vol:t:r:r}.WM.$leftprefix.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_binarize --i  $workdir/${vol:t:r:r}.WM.$leftprefix.mgz --o $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.mgz --min 1 --max 1 --binval $leftlabel)
  echo $cmd; eval $cmd

  # Right
  set filelist = ()
  foreach label (41 48 49 50 51 52 57 58 60 62) 
    set cmd = (mri_binarize --i $vol --o $workdir/${vol:t:r:r}.$rightprefix.$label.mgz --min $label --max $label --binval $rightlabel)
    foreach temp ($alllabels) 
      if ($temp == $label) then
        set filelist = ($filelist $workdir/${vol:t:r:r}.$rightprefix.$label.mgz)
        if !(-e $workdir/${vol:t:r:r}.$rightprefix.$label.mgz) then
          echo $cmd; eval $cmd
        endif
        break
      endif
    end
  end

  set cmd = (mri_or  $filelist $workdir/${vol:t:r:r}.WM.$rightprefix.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_binarize --i  $workdir/${vol:t:r:r}.WM.$rightprefix.mgz --o $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.mgz --min 1 --max 1 --binval $rightlabel)
  echo $cmd; eval $cmd

  # Brainstem # TODO allow 16 here
  foreach label (173 174 175 16) 
    set cmd = (mri_binarize --i $vol --o $workdir/${vol:t:r:r}.$label.mgz --min $label --max $label --binval 1)
    foreach temp ($alllabels) 
      if ($temp == $label) then
        if !(-e $workdir/${vol:t:r:r}.$label.mgz) then
          echo $cmd; eval $cmd
        endif
        break
      endif
    end
  end

endif #  getwmhemis

###
###  filling the holes and putting the two hemis together in one volume: FilledWM.mgz 
###

if ($getfilledwm) then
  # First fill the ventricles (L:4, R:43)
  set label = 4
  set cmd = (mri_binarize --i $vol --o $workdir/${vol:t:r:r}.$leftprefix.$label.mgz --min $label --max $label --binval $leftlabel)
  foreach temp ($alllabels)
    if ($temp == $label) then
      if !(-e $workdir/${vol:t:r:r}.$leftprefix.$label.mgz) then
        echo $cmd; eval $cmd
      endif
      set cmd = (mri_or $workdir/${vol:t:r:r}.$leftprefix.$label.mgz $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.mgz $workdir/${vol:t:r:r}.WM.$leftprefix.vfilled.mgz)
      if !(-e $workdir/${vol:t:r:r}.WM.$leftprefix.vfilled.mgz) then 
        echo $cmd; eval $cmd
      endif
      break
    endif
  end
  if !(-e $workdir/${vol:t:r:r}.WM.$leftprefix.vfilled.mgz) then
    set cmd = (mri_binarize --i $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.mgz --o  $workdir/${vol:t:r:r}.WM.$leftprefix.vfilled.mgz --min $leftlabel --max $leftlabel --binval 1)
    echo $cmd; eval $cmd
  endif
  set cmd = (mri_binarize --i $workdir/${vol:t:r:r}.WM.$leftprefix.vfilled.mgz --o $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.vfilled.mgz --min 1 --max 1 --binval $leftlabel)
  echo $cmd; eval $cmd

  set label = 43
  set cmd = (mri_binarize --i $vol --o $workdir/${vol:t:r:r}.$rightprefix.$label.mgz --min $label --max $label --binval $rightlabel)
  foreach temp ($alllabels)
    if ($temp == $label) then
      if !(-e $workdir/${vol:t:r:r}.$rightprefix.$label.mgz) then
        echo $cmd; eval $cmd
      endif
      set cmd = (mri_or $workdir/${vol:t:r:r}.$rightprefix.$label.mgz $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.mgz $workdir/${vol:t:r:r}.WM.$rightprefix.vfilled.mgz) 
      if !(-e $workdir/${vol:t:r:r}.WM.$rightprefix.vfilled.mgz) then
        echo $cmd; eval $cmd
      endif
      break
    endif
  end
  if !(-e $workdir/${vol:t:r:r}.WM.$rightprefix.vfilled.mgz) then
    set cmd = (mri_binarize --i $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.mgz --o  $workdir/${vol:t:r:r}.WM.$rightprefix.vfilled.mgz --min $rightlabel --max $rightlabel --binval 1)
    echo $cmd; eval $cmd
  endif
  set cmd = (mri_binarize --i $workdir/${vol:t:r:r}.WM.$rightprefix.vfilled.mgz --o $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.vfilled.mgz --min 1 --max 1 --binval $rightlabel)
  echo $cmd; eval $cmd

  # right 
  set input = $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.vfilled.mgz
  set output = ${input:r:r}.filled.mgz
  set cmd = (mri_morphology $input fill_holes 26 $output)
  echo $cmd; eval $cmd
  set cmd = (mri_extract_largest_CC -hemi rh -T 1 $output ${output:r}.CC.mgz)
  echo $cmd; eval $cmd

  # left
  set input = $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.vfilled.mgz
  set output = ${input:r:r}.filled.mgz
  set cmd = (mri_morphology $input fill_holes 26 $output)
  echo $cmd; eval $cmd
  set cmd = (mri_extract_largest_CC -hemi lh -T 1 $output ${output:r}.CC.mgz)
  echo $cmd; eval $cmd

  set cmd = (mri_or -o $workdir/${vol:t:r:r}.WM.$rightprefix.$rightlabel.filled.CC.mgz $workdir/${vol:t:r:r}.WM.$leftprefix.$leftlabel.filled.CC.mgz $workdir/FilledWM.mgz)
  echo $cmd; eval $cmd

  set cmd = ( cp $workdir/FilledWM.mgz $workdir/FilledWM.orig.mgz)
  echo $cmd; eval $cmd;
  set cmd = ( mri_morphology $workdir/FilledWM.mgz close 1 $workdir/FilledWM.mgz ) 
  echo $cmd; eval $cmd;

endif # getfilledwm

### 
### Getting ready for constructing surfaces  
### 

if ($surfaceprep) then
  # right
  mkdir -p $workdir/surf/

  set hemi = rh
  set cmd = (mri_pretess $workdir/FilledWM.mgz $rightlabel $mridir/norm.nii.gz $workdir/FilledWM-pretess$rightlabel.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_tessellate $workdir/FilledWM-pretess$rightlabel.mgz $rightlabel $workdir/surf/$hemi.orig.nofix)
  echo $cmd; eval $cmd

  set cmd = (mris_extract_main_component $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.orig.nofix)
  echo $cmd; eval $cmd

  set cmd = (mris_smooth -nw $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.smoothorig.nofix)
  echo $cmd; eval $cmd
  set cmd = (mris_inflate -no-save-sulc $workdir/surf/$hemi.smoothorig.nofix $workdir/surf/$hemi.inflated.nofix)
  echo $cmd; eval $cmd
  set cmd = (mris_sphere -q -in 3000 $workdir/surf/$hemi.inflated.nofix $workdir/surf/$hemi.qsphere.nofix) 
  echo $cmd; eval $cmd
  cp $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.orig  
  cp $workdir/surf/$hemi.inflated.nofix $workdir/surf/$hemi.inflated 

  # left
  set hemi = lh
  set cmd = (mri_pretess $workdir/FilledWM.mgz $leftlabel $mridir/norm.nii.gz $workdir/FilledWM-pretess$leftlabel.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_tessellate $workdir/FilledWM-pretess$leftlabel.mgz $leftlabel $workdir/surf/$hemi.orig.nofix)
  echo $cmd; eval $cmd

  set cmd = (mris_extract_main_component $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.orig.nofix)
  echo $cmd; eval $cmd 

  set cmd = (mris_smooth -nw $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.smoothorig.nofix)
  echo $cmd; eval $cmd
  set cmd = (mris_inflate -no-save-sulc $workdir/surf/$hemi.smoothorig.nofix $workdir/surf/$hemi.inflated.nofix)
  echo $cmd; eval $cmd
  set cmd = (mris_sphere -q -in 3000 $workdir/surf/$hemi.inflated.nofix $workdir/surf/$hemi.qsphere.nofix) 
  echo $cmd; eval $cmd
  cp $workdir/surf/$hemi.orig.nofix $workdir/surf/$hemi.orig  
  cp $workdir/surf/$hemi.inflated.nofix $workdir/surf/$hemi.inflated 

  pushd $surfdir 
  foreach hemi (lh rh)
    if (-e $hemi.qsphere.nofix) then 
      rm -f $hemi.qsphere.nofix
    endif
    set cmd = (cp $workdir/surf/$hemi.qsphere.nofix $hemi.qsphere.nofix)
    echo $cmd; eval $cmd
  end
  popd

  ###
  ### wm.mgz
  ###
  if (-e $workdir/${segfile:t:r:r}.173.mgz) then 
    set cmd = (mri_or $workdir/${segfile:t:r:r}.$rightprefix.41.mgz $workdir/${segfile:t:r:r}.$leftprefix.2.mgz \
                      $workdir/${segfile:t:r:r}.173.mgz $workdir/${segfile:t:r:r}.174.mgz \
                      $workdir/${segfile:t:r:r}.175.mgz $workdir/wm.nonscaled.mgz )
  else
    set cmd = (mri_or $workdir/${segfile:t:r:r}.$rightprefix.41.mgz $workdir/${segfile:t:r:r}.$leftprefix.2.mgz \
                      $workdir/${segfile:t:r:r}.16.mgz $workdir/wm.nonscaled.mgz )
  endif
  echo $cmd; eval $cmd
  set cmd = (mri_binarize --i $workdir/wm.nonscaled.mgz  --o $workdir/wm.mgz --min 1 --max 1 --binval 110)
  echo $cmd; eval $cmd

  #
  pushd $mridir
  if (-e wm.mgz) then 
    rm -f wm.mgz
  endif

  set cmd = (mri_convert --no_scale 1 --out_data_type uchar $workdir/wm.mgz $workdir/wm.uchar.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_edit_wm_with_aseg $workdir/wm.uchar.mgz brain.nii.gz $vol $workdir/wm.asegedit.mgz)
  echo $cmd; eval $cmd
  set cmd = (mri_pretess $workdir/wm.asegedit.mgz wm norm.nii.gz wm.mgz)
  echo $cmd; eval $cmd
  popd 
  #
  pushd $surfdir
  foreach hemi (lh rh)
    if (-e $hemi.orig) then 
      rm -f $hemi.orig
    endif
    set cmd = (cp $workdir/surf/$hemi.orig $hemi.orig)
    echo $cmd; eval $cmd
    if (-e $hemi.inflated) then 
      rm -f $hemi.inflated
    endif
    set cmd = (cp $workdir/surf/$hemi.inflated $hemi.inflated)
    echo $cmd; eval $cmd
    if (-e $hemi.qsphere) then 
      rm -f $hemi.qsphere
    endif
    set cmd = (cp $hemi.qsphere.nofix $hemi.qsphere)
    echo $cmd; eval $cmd

    set cmd = (mris_euler_number $hemi.orig) # 06/14/2017
    $cmd >& $hemi.orig.euler.txt
    set defect_idx = `grep "total defect index" $hemi.orig.euler.txt | awk '{print $5}'`
    if ($defect_idx > 0) then 
      set cmd = (mris_topo_fixer -mgz -warnings $subj $hemi)  # produces: $surfdir/$hemi.orig_corrected
      echo $cmd; eval $cmd
    else
      cp $hemi.orig $hemi.orig_corrected
    endif

    if (-e $surfdir/$hemi.orig) then 
      rm -f $surfdir/$hemi.orig
    endif
    set cmd = (cp $surfdir/$hemi.orig_corrected $surfdir/$hemi.orig)
    echo $cmd; eval $cmd
    set cmd = (mris_euler_number $surfdir/$hemi.orig )
    echo $cmd; eval $cmd
    set cmd = (mris_remove_intersection $surfdir/$hemi.orig $surfdir/$hemi.orig) 
    echo $cmd; eval $cmd
  end
  popd
endif # surfaceprep

###
if ($getwhitesurfaces) then

  pushd $mridir

  set cmd = (mri_mask -T 5 $mridir/brain.mgz $mridir/brainmask.mgz $mridir/brain.finalsurfs.mgz )
  echo $cmd; eval $cmd

  if (-e $mridir/filled.mgz) then 
    rm -f $mridir/filled.mgz
  endif
  #ln -s $workdir/FilledWM.mgz $mridir/filled.mgz
  cp $workdir/FilledWM.mgz $mridir/filled.mgz

  foreach hemi (lh rh)
    set cmd = (mris_make_surfaces -output .dist -soap -orig_white orig -aseg aseg \
               -cover_seg $vol -noaparc -whiteonly -mgz -T1 brain.finalsurfs $subj $hemi)
    echo $cmd; eval $cmd
    if (-e $surfdir/$hemi.white) then 
      rm -f $surfdir/$hemi.white
    endif
    set cmd = (mv $surfdir/$hemi.white.dist $surfdir/$hemi.white)
    echo $cmd; eval $cmd
    set cmd = (mv $surfdir/$hemi.curv.dist $surfdir/$hemi.curv)
    echo $cmd; eval $cmd
    set cmd = (mv $surfdir/$hemi.area.dist $surfdir/$hemi.area)
    echo $cmd; eval $cmd
    set cmd = (mv $labeldir/$hemi.cortex.dist.label $labeldir/$hemi.cortex.label)
    echo $cmd; eval $cmd
  end
  
  ###
  ### Smooting and other measures / representations
  ###
  foreach hemi (lh rh)
    
    set cmd = (mris_smooth -n 3 -nbrs 1 -gt .995 -d 2 -nw $surfdir/$hemi.white $surfdir/$hemi.smoothwm1) 
    echo $cmd; eval $cmd
    set cmd = (mris_smooth -n 3 -nbrs 1 -gt .995 -d 2 -nw $surfdir/$hemi.smoothwm1 $surfdir/$hemi.smoothwm2) 
    echo $cmd; eval $cmd
    set cmd = (mris_smooth -n 3 -nbrs 1 -gt .995 -d 2 -nw $surfdir/$hemi.smoothwm2 $surfdir/$hemi.smoothwm3) 
    echo $cmd; eval $cmd
    set cmd = (mris_smooth -n 3 -nbrs 1 -gt .995 -d 2 -nw $surfdir/$hemi.smoothwm3 $surfdir/$hemi.smoothwm4) 
    echo $cmd; eval $cmd
    set cmd = (mris_smooth -n 3 -nbrs 1 -gt .995 -d 2 -nw $surfdir/$hemi.smoothwm4 $surfdir/$hemi.smoothwm5) 
    echo $cmd; eval $cmd
    cp $surfdir/$hemi.smoothwm5 $surfdir/$hemi.white
    set cmd = (mris_smooth -n 3 -nw $surfdir/$hemi.white $surfdir/$hemi.smoothwm) 
    echo $cmd; eval $cmd

    # .sulc is also produced!
    set cmd = (mris_inflate $surfdir/$hemi.smoothwm $surfdir/$hemi.inflated) 
    echo $cmd; eval $cmd
    # inlated.H/.K
    set cmd = (mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 $surfdir/$hemi.inflated) 
    echo $cmd; eval $cmd 
    set cmd = (mris_sphere $surfdir/$hemi.inflated $surfdir/$hemi.sphere) 
    echo $cmd; eval $cmd

    # TODO: new steps 
    pushd $surfdir
    set cmd = (mris_register -curv $hemi.sphere $FREESURFER_HOME/average/$hemi.average.curvature.filled.buckner40.tif $hemi.sphere.reg)
    echo $cmd; eval $cmd

    set cmd = (mris_ca_label -l $labeldir/$hemi.cortex.label -aseg $mridir/aseg.mgz $subj $hemi $hemi.sphere.reg \
                             $FREESURFER_HOME/average/$hemi.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs $hemi.aparc.annot )
    echo $cmd; eval $cmd

    set cmd = (mris_ca_label -l $labeldir/$hemi.cortex.label -aseg $mridir/aseg.mgz $subj $hemi $hemi.sphere.reg \
                             $FREESURFER_HOME/average/$hemi.destrieux.simple.2009-07-29.gcs $hemi.aparc.a2009s.annot )
    echo $cmd; eval $cmd
    popd
  end
endif # getwhitesurfaces  

