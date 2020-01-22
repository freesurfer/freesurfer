#! /bin/tcsh -ef
#set echo=1

# this version of the original do,qbi script works with the Bay 3 SSFP diffusion data. 
# It includes the ability to reorient hemi scans, as well as whole brains, 
# alters the eddy parameters to work with SSFP diffusion data 
# and incorporates Benner's temperature drift correction MATLAB script

# 11/27/2017
# last edits include changes by lzollei in order to fix gradient table misorientation 
# issues
# ** some important updates: computing both distance normalized and unnormalized versions and probab. 
# curvature threshold differs between deterministic and probabilitic reconstructions!! (in order to match probab default) 


if ($#argv == 0) then
  echo "USAGE: $0 configfile"
  exit 1
endif

echo "Sourcing configuration from $argv[1]"
source $argv[1]
echo "properly sourced configuration file"
if (! $?dcmdir) then
  echo "ERROR: Must specify dcmdir in configuration file"
endif

if (! $?dcmlist) then
  echo "ERROR: Must specify dcmlist in configuration file"
endif

# if (! $?protocol) then
#  echo "ERROR: Must specify protocol in configuration file"
# endif

if (! $?protocol) then
  if (! $?pe1) then
    echo "ERROR: Must specify protocol or pe1 in configuration file"
  endif

  if (! $?pe2) then
    echo "ERROR: Must specify protocol or pe2 in configuration file"
  endif

  if (! $?rotime) then
    echo "ERROR: Must specify protocol or rotime in configuration file"
  endif
endif

if (! $?bvecfile) then
  echo "ERROR: Must specify bvecfile in configuration file"
endif

if (! $?bvalfile) then
  echo "ERROR: Must specify bvalfile in configuration file"
endif

if (! $?dwidir) then
  echo "ERROR: Must specify dwidir in configuration file"
endif

if (! $?doconcat) then	
  set doconcat = 0
endif
if (! $?dodrift) then	
  set dodrift = 0
endif
if (! $?doorient) then	
  set doorient = 0
endif
if (! $?domask) then	
  set domask = 0
endif
if (! $?doeddy) then	
  set doeddy = 0
endif
if (! $?dotensor) then	
  set dotensor = 0
endif
if (! $?doodf) then	
  set doodf = 0
endif
if (! $?dotrk) then 	
  set dotrk = 0
endif
if (! $?doprebed) then	
  set doprebed = 0
endif
if (! $?dobed) then	
  set dobed = 0
endif
if (! $?dopostbed) then	
  set dopostbed = 0
endif
if (! $?doprobtrk) then 
  set doprobtrk = 0
endif
if (! $?dopd) then
  set dopd = 0
endif

if (! $?lobmaskthresh) then	
  set lobmaskthresh = 0.1
endif
if (! $?hibmaskthresh) then	
  set hibmaskthresh = 400
endif  
if (! $?usehibmask) then	
  set usehibmask = 0
endif
if (! $?angthresh) then		
  set angthresh = 60   
endif
if (! $?nstick) then		
  set nstick = 2
endif
if (! $?seedroidir) then	
  set seedroidir = $dwidir/rois
endif 

if ($usehibmask) then
  set mask = highb
else
  set mask = lowb 
endif
 

echo "Everything from the configuration file has been properly sourced."
 

# Create output directory
 umask 002
 mkdir -p $dwidir
 set LF = $dwidir/log.txt


#
# Find b=0 and b>0 volumes
#
 set lowblist = ()
 set highblist = ()
 
set bvals = `cat $bvalfile`
# echo $bvals

 if (! $?keepframe) then
  set keepframe = `printf "%d\n" $bvals | awk '{print 1}'`
 endif
# echo $keepframe
 
 @ k = 1
 while ($k <= $#bvals)
  echo $k
  if ($bvals[$k] == 0) then
    if ($keepframe[$k]) set lowblist  = ($lowblist  `echo "$k-1" | bc`)
  else
    if ($keepframe[$k]) set highblist = ($highblist `echo "$k-1" | bc`)
   endif
   @ k = $k + 1
 end
 
if ($doconcat) then
  #
  # Convert sets of DWIs from dicom
  #
  set dwilist = ()

  @ k = 1
  while ($k <= $#dcmlist)
    set dcmfile = `grep trufi_diff $dcmdir/scan.log \
                   | awk -v run=$dcmlist[$k] '{if ($1 == run) print $NF}'`

    set dwifile = dwi_set`printf '%02d\n' $k`.nii.gz

    setenv FS_SAME_SLICE_THRESH .01 # needed in order to avoid strange exvivo save/trasfer error (LZ:08/09/2017)

    set cmd = mri_convert
    set cmd = ($cmd $dcmdir/$dcmfile)
    set cmd = ($cmd $dwidir/$dwifile)
    echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF

    set dwilist = ($dwilist $dwidir/$dwifile)

    @ k = $k + 1
  end

  #
  # Concatenate DWI sets
  #
  set cmd = mri_concat
  set cmd = ($cmd --i $dwilist)
  set cmd = ($cmd --o $dwidir/dwi_orig.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($dodrift) then
  #
  # Normalize images to compensate for temperature drift
  #
#  set cmd = "addpath $FREESURFER_HOME/matlab"
#  set cmd = "addpath /autofs/space/turan_001/users/lzollei/dev/matlab"
  set cmd = "addpath /usr/pubsw/packages/matlab/R2015b/bin/matlab"  # 01/22/2020: changed due to consistent failure when processing exvivo data -- not really sure what in the newer version makes things brake
  set cmd = "$cmd; fix_exvivo_dwi_drift("
  set cmd = "$cmd '$dwidir/dwi_drift.nii.gz', "
  set cmd = "$cmd '$dwidir/drift', "
  set cmd = "$cmd '$dwidir/dwi_orig.nii.gz', "
  set cmd = "$cmd  $#lowblist);"
  echo $cmd | tee -a $LF
  echo $cmd | matlab -nosplash -nodesktop
  # echo $cmd | /usr/pubsw/bin/matlab9.5 -nosplash
endif

#############
##### NOTE: given that the exvivo diffusion matrices are loaded from outside of the 
#####       images' dicom files they need to be moved into the image domain
#####       Here bvecs are rotated by the abs value of the ras2vox matrix
#####       This needs to happen before any of the operations with bvecs start
## 07/11/2017 (LZ)
set cmd = (mri_info --ras2vox $dwidir/dwi_orig.nii.gz)  
echo $cmd | tee -a $LF
eval $cmd  >> $dwidir/mtx.r2v.txt

if (-e $dwidir/neg.mtx.r2v.txt) then 
  rm -f $dwidir/neg.mtx.r2v.txt
endif 
foreach k (1 2 3 4)
  set line = `sed -n "$k p"  $dwidir/mtx.r2v.txt`
  echo `echo "$line[1]" | awk ' { if($1>=0) { print $1} else {print $1*-1 }}'` `echo "$line[2]" | awk ' { if($1>=0) { print $1} else {print $1*-1 }}'` `echo "$line[3]" | awk ' { if($1>=0) { print $1} else {print $1*-1 }}'` `echo "$line[4]" | awk ' { if($1>=0) { print $1} else {print $1*-1 }}'` >> $dwidir/neg.mtx.r2v.txt
end
# more $dwidir/neg.mtx.r2v.txt

set cmd = (xfmrot $dwidir/neg.mtx.r2v.txt $bvecfile $dwidir/r2v.bvecs)
echo $cmd | tee -a $LF
$cmd |& tee -a $LF
set bvecfile = $dwidir/r2v.bvecs
#############

if ($doorient) then
  if (-e $dwidir/dwi_drift.nii.gz) then
    set dwiname = dwi_drift
  else
    set dwiname = dwi_orig
  endif

  #
  # Rotate by 90 degrees around the x-axis to get the true brain orientation:
  # (x, y, z) -> (x, -z, y)
  # Assuming that the ex vivo brain is placed with the brainstem down
  # and the frontal lobe pointing towards the console room,
  # the axis that one would expect to be S-I in a live human
  # is instead occupied by the P-A axis of the ex vivo brain 
  # 
  # NOTE: If the brain were placed with the occipital lobe pointing
  # towards the console room, an additional rotation around the y-axis
  # would be needed that is not performed here!
  #

  set cdc =  `mri_info --cdc $dwidir/$dwiname.nii.gz | awk '{print $1, -$3, $2}'`
  set rdc =  `mri_info --rdc $dwidir/$dwiname.nii.gz | awk '{print $1, -$3, $2}'`
  set sdc =  `mri_info --sdc $dwidir/$dwiname.nii.gz | awk '{print $1, -$3, $2}'`
  set cras = `mri_info --cras $dwidir/$dwiname.nii.gz | awk '{print $1, -$3, $2}'`

  if ($?hemi) then
    if ($hemi == lh) then
      #
      # This sample is a left hemisphere, placed lateral side down
      # Rotate by 90 degrees around the y-axis to get the true brain orientation:
      # (x, y, z) -> (z, y, -x)
      #
      set cdc =  `echo $cdc | awk '{print $3, $2, -$1}'`
      set rdc =  `echo $rdc | awk '{print $3, $2, -$1}'`
      set sdc =  `echo $sdc | awk '{print $3, $2, -$1}'`
      set cras = `echo $cras | awk '{print $3, $2, -$1}'`

    else if ($hemi == rh) then
      #
      # This sample is a right hemisphere, placed lateral side down
      # Rotate by -90 degrees around the y-axis to get the true brain orientation:
      # (x, y, z) -> (-z, y, x)
      #
      set cdc =  `echo $cdc | awk '{print -$3, $2, $1}'`
      set rdc =  `echo $rdc | awk '{print -$3, $2, $1}'`
      set sdc =  `echo $sdc | awk '{print -$3, $2, $1}'`
      set cras = `echo $cras | awk '{print -$3, $2, $1}'`

    endif
  endif


  set cmd = mri_convert
  set cmd = ($cmd -iid $cdc)
  set cmd = ($cmd -ijd $rdc)
  set cmd = ($cmd -ikd $sdc)
  set cmd = ($cmd -ic  $cras)
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd $dwidir/$dwiname.hdrorient.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (cp $bvecfile $dwidir/$dwiname.hdrorient.bvecs)
  echo $cmd; eval $cmd
  echo $cmd |& tee -a $LF

  set cmd = (cp $bvalfile $dwidir/$dwiname.hdrorient.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  #  Convert DWIs and gradient vectors to LAS orientation
  #  (Use this instead of flip4fsl, which inverts bvecs in x erroneously)
  #  Note, this will also produce the reoriented bvec file 
  #  NB: bvec file needs to be consistently named with corresponding dwi!!!
  set cmd = orientLAS
  set cmd = ($cmd $dwidir/$dwiname.hdrorient.nii.gz)
  set cmd = ($cmd $dwidir/dwi_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif


if ($domask) then
  #
  # Extract brain mask from low-b images
  #
  set cmd = mri_convert
  set cmd = ($cmd $dwidir/dwi_las.nii.gz)
  set cmd = ($cmd -f $lowblist)
  set cmd = ($cmd $dwidir/lowb_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = mri_concat
  set cmd = ($cmd --i $dwidir/lowb_las.nii.gz)
  set cmd = ($cmd --mean)
  set cmd = ($cmd --o $dwidir/lowb_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF


  set cmd = bet
  set cmd = ($cmd $dwidir/lowb_las.nii.gz)
  set cmd = ($cmd $dwidir/lowb_las_brain.nii.gz)
  set cmd = ($cmd -m -f $lobmaskthresh)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  # Extract mask from high-b images
  # This is done in case mask from low-b fails because the brain is in formalin
  # instead of fomblin, so there is a lot of background
  #
  set cmd = mri_convert
  set cmd = ($cmd $dwidir/dwi_las.nii.gz)
  set cmd = ($cmd -f $highblist)
  set cmd = ($cmd $dwidir/highb_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = mri_concat
  set cmd = ($cmd --i $dwidir/highb_las.nii.gz)
  set cmd = ($cmd --mean)
  set cmd = ($cmd --o $dwidir/highb_las.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = mri_binarize
  set cmd = ($cmd --i $dwidir/highb_las.nii.gz)
  set cmd = ($cmd --min $hibmaskthresh --dilate 4 --erode 3)
  set cmd = ($cmd --o $dwidir/highb_las_brain_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($doeddy) then
  #
  # Compensate for eddy-current distortions
  #
  if ($?protocol) then	# Read acquisition parameters from protocol print-out
 			# echo spacing and EPI factor not relevant to SSFP scan
    set pe1 = `awk '{if ($1 ~ /trufi_diff_tb/) \
                        {while ($0 !~ /Phase enc. dir./) {getline}; pe1 = $4}} \
                    END {print pe1}' $protocol`

    set pe2 = `awk '{if ($1 ~ /trufi_diff_tb/) \
                        {while ($0 !~ /Phase enc. dir./) {getline}; pe2 = $6}} \
                    END {print pe2}' $protocol`

    #set esp = `awk '{if ($1 ~ /trufi_diff_tb/) \
    #                    {while ($0 !~ /Echo spacing/) {getline}; esp = $3}} \
    #                END {print esp}' $protocol`

    #set epif = `awk '{if ($1 ~ /trufi_diff_tb/) \
    #                    {while ($0 !~ /EPI factor/) {getline}; epif = $3}} \
    #                 END {print epif}' $protocol`

    set rotime = 0.1 # 1 doesn't work in eddy
    #    set rotime = `echo "$esp * 0.001 * ($epif - 1)" | bc -l`
  endif

  if ($pe1 == R && $pe2 == L) then
    echo "1 0 0 $rotime"    > $dwidir/acqp.txt
  else if ($pe1 == L && $pe2 == R) then
    echo "-1 0 0 $rotime"   > $dwidir/acqp.txt
  else if ($pe1 == P && $pe2 == A) then
    echo "0 1 0 $rotime"    > $dwidir/acqp.txt
  else if ($pe1 == A && $pe2 == P) then
    echo "0 -1 0 $rotime"   > $dwidir/acqp.txt
  endif

  set bvals = `cat $dwidir/dwi_las.bvals`
  echo `printf "%d\n" $bvals | awk '{print 1}'` > $dwidir/index.txt

  set cmd = eddy
  set cmd = ($cmd --imain=$dwidir/dwi_las.nii.gz)
  set cmd = ($cmd --mask=$dwidir/${mask}_las_brain_mask.nii.gz)
  set cmd = ($cmd --acqp=$dwidir/acqp.txt)
  set cmd = ($cmd --index=$dwidir/index.txt)
  set cmd = ($cmd --bvecs=$dwidir/dwi_las.bvecs)
  set cmd = ($cmd --bvals=$dwidir/dwi_las.bvals)
  set cmd = ($cmd --out=$dwidir/dwi)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  # Apply rotations from eddy-current compensation to gradient vectors
  #
  set cmd = xfmrot
  set cmd = ($cmd $dwidir/dwi.eddy_parameters)
  set cmd = ($cmd $dwidir/dwi_las.bvecs)
  set cmd = ($cmd $dwidir/dwi.bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (cp $dwidir/dwi_las.bvals $dwidir/dwi.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  # Extract brain mask from average low-b image (corrected)
  #
  set cmd = (mri_convert $dwidir/dwi.nii.gz -f $lowblist $dwidir/lowb.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (mri_concat --i $dwidir/lowb.nii.gz --mean --o $dwidir/lowb.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = bet
  set cmd = ($cmd $dwidir/lowb.nii.gz)
  set cmd = ($cmd $dwidir/lowb_brain.nii.gz)
  set cmd = ($cmd -m -f $lobmaskthresh)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  #
  # Extract mask from high-b images (corrected)
  # This is done in case mask from low-b fails because the brain is in formalin
  # instead of fomblin, so there is a lot of background
  #
  set cmd = (mri_convert $dwidir/dwi.nii.gz -f $highblist $dwidir/highb.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = mri_concat
  set cmd = ($cmd --i $dwidir/highb.nii.gz)
  set cmd = ($cmd --mean)
  set cmd = ($cmd --o $dwidir/highb.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = mri_binarize
  set cmd = ($cmd --i $dwidir/highb.nii.gz)
  set cmd = ($cmd --min $hibmaskthresh --dilate 4 --erode 3)
  set cmd = ($cmd --o $dwidir/highb_brain_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif


if ($dotensor) then
  #
  # Fit tensors
  #
  if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
    set maskname = $mask
  else
    set dwiname = dwi_las
    set maskname = ${mask}_las
  endif

  set cmd = dtifit
  set cmd = ($cmd -k $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -o $dwidir/dtifit)
  set cmd = ($cmd -m $dwidir/${maskname}_brain_mask.nii.gz)
  set cmd = ($cmd -r $dwidir/$dwiname.bvecs)
  set cmd = ($cmd -b $dwidir/$dwiname.bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

set dtkdir = /usr/pubsw/packages/dtk/0.6.4.1_patched

if ($doodf) then
  #
  # Fit ODFs
  #
  if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
  else
    set dwiname = dwi_las
  endif


  # Move low-b images to the beginning of the series (for dtk)
  set cmd = mri_convert
  set cmd = ($cmd $dwidir/$dwiname.nii.gz)
  set cmd = ($cmd -f $lowblist $highblist)
  set cmd = ($cmd $dwidir/tmp.$dwiname.nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  cp /dev/null $dwidir/gradients.txt


  foreach k ($highblist)
    set k1 = `echo "$k+1" | bc`
    set xyz = `sed -n "$k1 p" $dwidir/$dwiname.bvecs`
    echo "$xyz[1], $xyz[2], $xyz[3]" >> $dwidir/gradients.txt
  end

  set cmd = $dtkdir/hardi_mat
  set cmd = ($cmd $dwidir/gradients.txt)
  set cmd = ($cmd $dwidir/qbi_mat.dat)
  set cmd = ($cmd -ref $dwidir/tmp.$dwiname.nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set nlow = $#lowblist
  set ndir = `echo "$#highblist + 1" | bc`

  setenv DSI_PATH $dtkdir/matrices

  set cmd = $dtkdir/odf_recon
  set cmd = ($cmd $dwidir/tmp.$dwiname.nii)
  set cmd = ($cmd $ndir 181)
  set cmd = ($cmd $dwidir/qbi)
  set cmd = ($cmd -b0 $nlow)
  set cmd = ($cmd -mat $dwidir/qbi_mat.dat)
  set cmd = ($cmd -nt -p 3 -sn 1 -ot nii)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
  
  #
  # Run dti tractography (Brian doesn't use dti track file, but does use FA, ADC and color FA maps generated by this)
  #

  set cmd = $dtkdir/dti_recon
  set cmd = ($cmd $dwidir/dwi.nii.gz)
  set cmd = ($cmd $dwidir/dti)
  set cmd = ($cmd -gm $dwidir/$dwiname.bvecs)
  set cmd = ($cmd -b 4080)
  set cmd = ($cmd -b0 $nlow)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF 

  # Clean up temporary files
  set cmd = (rm -f $dwidir/tmp.$dwiname.nii $dwidir/qbi_mat.dat)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($dotrk) then
  #
  # The foreach loop below run various inversions and swaps so you can check which version is correct
  # Run q-ball tractography
  #

  setenv DSI_PATH $dtkdir/matrices

  ##  Brian said no inversion, and yz swap is correct one, so following line should be used:
  ##  set trkfile = qbi.inv,no.swap,yz.trk  
  #
  # LZ (11/27/2017): the above is outdated - it is either (swap no, inv z) that works or (swap no, inv no) with a newer / correct version of DTK -- hard to track when that version got released and what version was used...
  #

  foreach inv (no x y z)
   foreach swap (no xy yz zx)
     set trkfile = qbi.inv,$inv.swap,$swap.trk
     set cmd = $dtkdir/odf_tracker
     set cmd = ($cmd $dwidir/qbi)
     set cmd = ($cmd $dwidir/tmp.qbi.trk)
     set cmd = ($cmd -at $angthresh)
     set cmd = ($cmd -m $dwidir/qbi_dwi.nii)
     set cmd = ($cmd -it nii)
     if ($inv != no)	set cmd = ($cmd -i$inv)
     if ($swap != no)	set cmd = ($cmd -s$swap)
     echo $cmd |& tee -a $LF
     $cmd |& tee -a $LF
	    
     set cmd = ($dtkdir/spline_filter $dwidir/tmp.qbi.trk 1 $dwidir/$trkfile)
     echo $cmd |& tee -a $LF
     $cmd |& tee -a $LF
	    
    # Clean up temporary files
     set cmd = (rm -f $dwidir/tmp.qbi.trk)
     echo $cmd |& tee -a $LF
    $cmd |& tee -a $LF
    end
  end
endif


if ($doprebed) then
  #
  # Do bedpostx pre-processing
  # This is done separately because it needs to be run on a machine
  # with more memory than launchpad nodes
  #
  if (-e $dwidir/dwi.nii.gz) then
    set dwiname = dwi
    set maskname = ${mask}
  else
    set dwiname = dwi_las
    set maskname = ${mask}_las
  endif

  set cmd = (ln -sf $dwidir/$dwiname.nii.gz $dwidir/data.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/$dwiname.bvals $dwidir/bvals)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/$dwiname.bvecs $dwidir/bvecs)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  set cmd = (ln -sf $dwidir/${maskname}_brain_mask.nii.gz)
  set cmd = ($cmd $dwidir/nodif_brain_mask.nii.gz)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF

  mkdir -p $dwidir.bedpostX
  mkdir -p $dwidir.bedpostX/diff_slices
  mkdir -p $dwidir.bedpostX/logs
  mkdir -p $dwidir.bedpostX/xfms

  set cmd = (bedpostx_preproc.sh $dwidir)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($dobed) then
  #
  # Run bedpostx main processing (slice by slice)
  #
  set nslice = `mri_info --nslices $dwidir/lowb.nii.gz`
  set cmdfile = $dwidir.bedpostX/commands.txt
  rm -f $cmdfile

  @ k = 0
  while ($k < $nslice)
    set cmd = bedpostx_single_slice.sh
    set cmd = ($cmd $dwidir $k)
    set cmd = ($cmd --nf=$nstick --fudge=1 --bi=1000 --nj=1250)
    set cmd = ($cmd --se=25 --model=1 --cnonlinear)
    echo $cmd >> $cmdfile
    @ k = $k + 1
  end

  set cmd = fsl_sub_mgh
  set cmd = ($cmd -l $dwidir.bedpostX/logs)
  set cmd = ($cmd -m a -N bedpostx)
  set cmd = ($cmd -t $cmdfile)
  set cmd = ($cmd -a ppn=1,vmem=32gb)
  echo "Run on launchpad:" |& tee -a $LF
  echo $cmd |& tee -a $LF
endif

if ($dopostbed) then
  #
  # Do bedpostx post-processing
  # This is done separately because it needs to be run on a machine
  # with more memory than launchpad nodes
  #
  # NOTE, cleanup for prebed happens here!
  #
  set cmd = (bedpostx_postproc.sh $dwidir)
  echo $cmd |& tee -a $LF
  $cmd |& tee -a $LF
endif

if ($doprobtrk) then
  #
  # Probabilistic seed-based targeted tractography -- all seeds are used as targets as well (LZ: 11/27/2017)
  #
  set seedlist = ($seedroidir/*.nii)

  echo "Found $#seedlist seed ROIs in $seedroidir"

  set step = `mri_info --cres $seedlist[1]`
  set step = `echo "$step / 3" | bc -l`
  set step = `printf '%g' $step`

  foreach seed ($seedlist)
    set seedname = `basename $seed | sed 's/.nii//'`

    foreach target ($seedlist)
      set targetname = `basename $target | sed 's/.nii//'`
      if !($seedname == $targetname) then
        set outdir = $dwidir/dpath.uncorrected.targeted.$seedname.2.$targetname # -pd is not on
        mkdir -p $outdir
      
        if !(-e $dwidir/rois/target$targetname.txt) then
          echo $dwidir/rois/$targetname.nii >> $dwidir/rois/target$targetname.txt
        endif 
 
        set cmd = probtrackx2
        set cmd = ($cmd -x $seed)
        set cmd = ($cmd -s $dwidir.bedpostX/merged)
        set cmd = ($cmd -m $dwidir.bedpostX/nodif_brain_mask)
        set cmd = ($cmd --dir=$outdir)
        set cmd = ($cmd -l --onewaycondition)
        set cmd = ($cmd -c 0.2 -S 2000 -P 5000) # Note: use -c .5 curvature threshold (vs .2, the default) in order to have correspondance with the 60 degrees angle threshold of deterministic tractography
        set cmd = ($cmd --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0)
        set cmd = ($cmd --steplength=$step)
        set cmd = ($cmd --forcedir --opd) 
        set cmd = ($cmd --targetmasks=$dwidir/rois/target$targetname.txt --stop=$dwidir/rois/target$targetname.txt --waypoints=$dwidir/rois/target$targetname.txt --os2t --s2tastext)
        echo $cmd |& tee -a $LF
        # $cmd |& tee -a $LF
        pbsubmit -m `whoami` -q p30 -c "$cmd" -l nodes=1:ppn=4,vmem=28gb



        set outdir = $dwidir/dpath.corrected.targeted.$seedname.2.$targetname # -pd is off while dopd=0 ("Correct path distribution for the length of the pathways")
        mkdir -p $outdir
	if ($dopd) then
  	  set cmd = ($cmd --pd)
	endif
        echo $cmd |& tee -a $LF
        # $cmd |& tee -a $LF
        pbsubmit -m `whoami` -q p30 -c "$cmd" -l nodes=1:ppn=4,vmem=28gb
      endif
    end
  end
endif
