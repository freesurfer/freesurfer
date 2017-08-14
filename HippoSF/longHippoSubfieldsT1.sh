#! /bin/tcsh -f

# If no arguments given
if( $#argv < 1 ) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   longHippoSubfieldsT1.sh BASE_SUBJECT_ID"
  echo " "
  echo "Or, for help"
  echo " "
  echo "   longHippoSubfieldsT1.sh --help"
  echo " "
  exit 1
endif

checkMCR
if($status) exit 1;

if( $1 == "--help") then
  echo " "
  echo "LONGITUDINAL HIPPOCAMPAL SUBFIELDS"
  echo " "
  echo "FreeSurfer currently supports longitudinal segmentation of the hippocampal"
  echo "subfields on T1-weighted images. The only requirement is that the time points"
  echo "do not mix 1mm scans processed with the standard stream and higher resolution"
  echo "scans processed with the -cm flag (conform to min voxel size for hi-res support)."
  echo " "
  echo "To use this module, you first need to run the FreeSurfer longitudinal stream"
  echo "on your data, i.e., steps 1-3 in Section LONGITUDINAL PROCESSING above: cross-"
  echo "sectional processing; creation and processing of base; and longitudinal processing."
  echo " "
  echo "Then, you can run the following command to segment the subfields from all time"
  echo "points simultaneously using the method described in [1]:"
  echo " "
  echo "   longHippoSubfieldsT1.sh BASE_SUBJECT_ID [SUBJECT_DIR]"
  echo " "
  echo "   (the argument [SUBJECT_DIR] is only necessary if the"
  echo "    environment variable has not been set)"
  echo " "
  echo "Note that you do not need to specify the IDs of the time points; those are read"
  echo "directly from the base subject directory"
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "http://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalHippocampalSubfields"
  echo " "
  echo "[1] Iglesias J.E., Van Leemput K., Augustinack J., Insausti R., Fischl B., Reuter M.,"
  echo "Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using"
  echo "subject-specific atlases, Neuroimage, in press."
  echo "http://dx.doi.org/10.1016/j.neuroimage.2016.07.020"
  echo " "
  exit 0
endif

# Error if SUBJECTS_DIR (the environment variable) does not exist
if ($#argv == 1) then
if (! $?SUBJECTS_DIR)  then
  echo " "
  echo "SUBJECTS_DIR variable does not exist"
  echo "Please define it or provide subjects directory as second input"
  echo " "
  exit 1
endif
endif

# Error if SUBJECTS_DIR (the environemnt variable) is empty 
if ($#argv == 1) then
if ( $SUBJECTS_DIR == "" ) then
  echo " "
  echo "SUBJECTS_DIR variable is empty"
  echo "Please redefine it or provide subjects directory as second input"
  echo " "
  exit 1
endif
endif

# If SUBJECTS_DIR is provided, just set it
if ($#argv == 2) then
  set SUBJECTS_DIR = $2
endif

# Set base subject
set BASESUBJ = $1

# Error if subject directory does not exist
if (! -d $SUBJECTS_DIR ) then
  echo " "
  echo "Subjects directory:"
  echo "   $SUBJECTS_DIR"
  echo "does not exist"
  echo " "
  exit 1
endif

# Error if base subject directory does not exist
if (! -d $SUBJECTS_DIR/$BASESUBJ ) then
  echo " "
  echo "Directory of base subject:"
  echo "   $SUBJECTS_DIR/$BASESUBJ"
  echo "does not exist"
  echo " "
  exit 1
endif

# Error if file with time points not found
if (! -e $SUBJECTS_DIR/$BASESUBJ/base-tps ) then
  echo " "
  echo "List of time points:"
  echo "   $SUBJECTS_DIR/$BASESUBJ/base-tps"
  echo "does not exist"
  echo " "
  exit 1
endif

# Error if base subject not processed far enough
if (! -e ${SUBJECTS_DIR}/${BASESUBJ}/mri/wmparc.mgz || \
    ! -e ${SUBJECTS_DIR}/${BASESUBJ}/mri/norm.mgz || \
    ! -e ${SUBJECTS_DIR}/${BASESUBJ}/mri/transforms/talairach.xfm ) then
  echo " "
  echo "Cannot find wmparc.mgz or norm.mgz or talairach.xfm for the base subject."
  echo " "
  exit 1;
endif 

# Go ahead and read time points
set SubjsList = (`cat $SUBJECTS_DIR/$BASESUBJ/base-tps`)

# Make sure that necessary files are in there
foreach s (${SubjsList})
  set longSubj = ${s}.long.$BASESUBJ
  
  if (! -e ${SUBJECTS_DIR}/${longSubj}/mri/wmparc.mgz || \
      ! -e ${SUBJECTS_DIR}/${longSubj}/mri/norm.mgz || \
      ! -e ${SUBJECTS_DIR}/${longSubj}/mri/transforms/talairach.xfm ) then
    echo "Cannot find wmparc.mgz or norm.mgz or talairach.xfm for $longSubj"
    exit 1;
  endif 
end


# Make sure that the longitudinal subfields are not running already for this subject
set IsRunningFile = ${SUBJECTS_DIR}/${BASESUBJ}/scripts/IsRunningLongHPsub.lh+rh
if(-e $IsRunningFile) then
    echo ""
    echo "It appears that longitudinal subfields is already running"
    echo "for this subject based on the presence of IsRunningLongHPsub.lh+rh."
    echo "It could also be that longitudinal subfields was running at one point"
    echo "but died in an unexpected way. If it is the case that there"
    echo "is a process running, you can kill it and start over or"
    echo "just let it run. If the process has died, you should type:"
    echo ""
    echo "rm $IsRunningFile"
    echo ""
    echo "and re-run. Or you can add -no-isrunning to the recon-all"
    echo "command-line. The contents of this file are:"
    echo "----------------------------------------------------------"
    cat  $IsRunningFile
    echo "----------------------------------------------------------"
    exit 1;
endif


# If everything is in place, let's do it! First, we create the IsRunning file
echo "------------------------------" > $IsRunningFile
echo "BASE SUBJECT $BASESUBJ" >> $IsRunningFile
echo "DATE `date`"     >> $IsRunningFile
echo "USER $user"      >> $IsRunningFile
echo "HOST `hostname`" >> $IsRunningFile
echo "PROCESSID $$ "   >> $IsRunningFile
echo "PROCESSOR `uname -m`" >> $IsRunningFile
echo "OS `uname -s`"       >> $IsRunningFile
uname -a         >> $IsRunningFile
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $IsRunningFile
endif


# Now the real job
set hippohemilist=(left right)
set HSFLOG = (${SUBJECTS_DIR}/${BASESUBJ}/scripts/long-hippocampal-subfields-T1.log)
rm -f $HSFLOG

echo "------------------------------" > $HSFLOG
echo "USER $user"      >> $HSFLOG
echo "HOST `hostname`" >> $HSFLOG
echo "PROCESSID $$ "   >> $HSFLOG
echo "PROCESSOR `uname -m`" >> $HSFLOG
echo "OS `uname -s`"       >> $HSFLOG
uname -a         >> $HSFLOG
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $HSFLOG
endif
echo "------------------------------" >> $HSFLOG
echo " " >> $HSFLOG

foreach hemi ($hippohemilist)
  echo "#--------------------------------------------" \
    |& tee -a $HSFLOG
  echo "#@# Longitudinal Hippocampal Subfields processing (T1) $hemi `date`" \
    |& tee -a $HSFLOG

  set cmd = "${FREESURFER_HOME}/bin/segmentSF_T1_long.sh  ${FREESURFER_HOME}/MCRv80 $FREESURFER_HOME $SUBJECTS_DIR $hemi $BASESUBJ"

  foreach s ($SubjsList)
    set cmd="$cmd ${s}.long.$BASESUBJ"
  end

  fs_time ls >& /dev/null
  if ($status) then
    $cmd |& tee -a $HSFLOG
  else
    fs_time $cmd |& tee -a $HSFLOG
  endif

  if ($status) then

    uname -a | tee -a $HSFLOG
    echo "" |& tee -a $HSFLOG
    echo "longitudinal subfields exited with ERRORS at `date`" \
    |& tee -a $HSFLOG
    echo "" |& tee -a $HSFLOG
    echo "For more details, see the log file $HSFLOG"
    echo ""
    exit 1;
  endif
end
 
# All done!
rm -f $IsRunningFile

echo ""
echo "If you have used results from this software for a publication, please cite:"
echo ""
echo "Iglesias J.E., Van Leemput K., Augustinack J., Insausti R., Fischl B., Reuter M.,"
echo "Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using"
echo "subject-specific atlases, Neuroimage, in press."
echo "http://dx.doi.org/10.1016/j.neuroimage.2016.07.020"
echo ""

exit 0



