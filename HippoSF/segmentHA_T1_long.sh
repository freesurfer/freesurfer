#! /bin/tcsh -f

set tcsh61706 = (`tcsh --version | grep "6\.17\.06"`)
if ("$tcsh61706" != "") then
  echo ""
  echo "WARNING: tcsh v6.17.06 has an exit code bug! Please update tcsh!"
  echo ""
  # workaround to force expected behavior:
  set anyerror
endif

# check MCR installation
checkMCR.sh
if ($status) then
    exit 1
endif

# If no arguments given
if( $#argv < 1 ) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   segmentHA_T1_long.sh BASE_SUBJECT_ID [SUBJECT_DIR]" 
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentHA_T1_long.sh --help"  
  echo " "
  exit 1
endif


if( $1 == "--help") then
  echo " "
  echo "LONGITUDINAL SEGMENTATION OF HIPPOCAMPAL SUBFIELDS AND NUCLEI OF THE AMYGDALA"
  echo " "
  echo "FreeSurfer currently supports longitudinal segmentation of the hippocampal"
  echo "subfields and nuclei of the amygdala on T1-weighted images. The only requirement"
  echo "is that the time points do not mix 1mm scans processed with the standard stream"
  echo "and higher resolution scans processed with the -cm flag"
  echo " "
  echo "To use this module, you first need to run the FreeSurfer longitudinal stream"
  echo "on your data, i.e., steps 1-3 in : "
  echo "https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing"
  echo "i.e., cross-sectional processing; creation and processing of base; and longitudinal"
  echo "processing."
  echo " "
  echo "Then, you can run the following command to segment the subfields/nuclei from all time"
  echo "points simultaneously using the method described in [1] (and the atlases described"
  echo "in [2] and [3]):"
  echo " "
  echo "   segmentHA_T1_long.sh BASE_SUBJECT_ID [SUBJECT_DIR]"
  echo " "
  echo "   (the argument [SUBJECT_DIR] is only necessary if the"
  echo "    environment variable SUBJECTS_DIR has not been set"
  echo "    or if you want to override it)"
  echo " "
  echo "Note that you do not need to specify the IDs of the time points; those are read"
  echo "directly from the base subject directory"
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "https://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalSubfieldsAndNucleiOfAmygdala"
  echo " "
  echo "[1] Iglesias J.E., Van Leemput K., Augustinack J., Insausti R., Fischl B., Reuter M.,"
  echo "Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using"
  echo "subject-specific atlases, Neuroimage, in press."
  echo "http://dx.doi.org/10.1016/j.neuroimage.2016.07.020"
  echo " "
  echo "[2] Iglesias, J.E., Augustinack, J.C., Nguyen, K., Player, C.M., Player, A., Wright,"
  echo "M., Roy, N., Frosch, M.P., McKee, A.C., Wald, L.L., Fischl, B., and Van Leemput, K.,"
  echo "A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution"
  echo "MRI: Application to adaptive segmentation of in vivo MRI.  Neuroimage 115, 2015, 117-137." 
  echo "http://dx.doi.org/10.1016/j.neuroimage.2015.04.042"
  echo " "
  echo "[3] Saygin, Z.M. & Kliemann, D. (joint 1st authors), Iglesias, J.E., van der Kouwe, A.J.W.,"
  echo "Boyd, E., Reuter, M., Stevens, A., Van Leemput, K., McKee, A., Frosch, M.P., Fischl, B.,"
  echo "and Augustinack, J.C., High-resolution magnetic resonance imaging reveals nuclei of the" 
  echo "human amygdala: manual segmentation to automatic atlas. Neuroimage 155, 2017, 370-382."
  echo "http://doi.org/10.1016/j.neuroimage.2017.04.046"
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
  set SUBJECTS_DIR = `getfullpath  $2`
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
  echo "Please make sure that the base subject has been processed."
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
    echo "Please make sure you have run the longitudinal processing for all subjects."
    exit 1;
  endif 
end


# Make sure that the longitudinal subfields are not running already for this subject
set IsRunningFile = ${SUBJECTS_DIR}/${BASESUBJ}/scripts/IsRunningLongHPsub.lh+rh
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that Longitudinal Subfields is already running"
  echo "for this subject based on the presence of IsRunningLongHPsub.lh+rh."
  echo "It could also be that Longitudinal Subfields was running at one point"
  echo "but died in an unexpected way. If it is the case that there"
  echo "is a process running, you can kill it and start over or"
  echo "just let it run. If the process has died, you should type:"
  echo ""
  echo "rm $IsRunningFile"
  echo ""
  echo "and re-run."
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

# Parameters
set RUNTIME="$FREESURFER_HOME/MCRv97/";
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K1="0.05";
set K2="0.05";
set OPTIMIZER="L-BFGS";
set MRFCONSTANT="0";
set SUFFIX="long.v22";

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

  # command
  set cmd="run_SegmentSubfieldsT1Longitudinal.sh $RUNTIME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K1 $K2 $hemi $OPTIMIZER $SUFFIX  '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/'  $MRFCONSTANT $BASESUBJ"

  foreach s ($SubjsList)
    set cmd="$cmd ${s}.long.$BASESUBJ"
  end

  fs_time ls >& /dev/null
  if ($status) then
    $cmd |& tee -a $HSFLOG
    set returnVal=$status
  else
    fs_time $cmd |& tee -a $HSFLOG
    set returnVal=$status
  endif

  if ($returnVal) then

    uname -a | tee -a $HSFLOG
    echo "" |& tee -a $HSFLOG
    echo "longitudinal subfields exited with ERRORS at `date`" \
    |& tee -a $HSFLOG
    echo "" |& tee -a $HSFLOG
    echo "For more details, see the log file $HSFLOG"
    echo ""
    rm -f $IsRunningFile
    exit 1;
  endif
end
 
# All done!
rm -f $IsRunningFile

echo " "
echo "All done!"
echo " "
echo "If you have used results from this software for a publication, please cite:"
echo ""
echo "Iglesias J.E., Van Leemput K., Augustinack J., Insausti R., Fischl B., Reuter M.,"
echo "Bayesian longitudinal segmentation of hippocampal substructures in brain MRI using"
echo "subject-specific atlases, Neuroimage, in press."
echo "http://dx.doi.org/10.1016/j.neuroimage.2016.07.020"
echo " "
echo "In addition, please cite the following paper if you have used the hippocampal subfields:"
echo " "
echo "Iglesias, J.E., Augustinack, J.C., Nguyen, K., Player, C.M., Player, A., Wright,"
echo "M., Roy, N., Frosch, M.P., McKee, A.C., Wald, L.L., Fischl, B., and Van Leemput, K.,"
echo "A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution"
echo "MRI: Application to adaptive segmentation of in vivo MRI.  Neuroimage 115, 2015, 117-137." 
echo "http://dx.doi.org/10.1016/j.neuroimage.2015.04.042"
echo " "
echo "And/or the following paper if you have used the segmentation of the nuclei of the amygdala:"
echo ""
echo "Saygin, Z.M. & Kliemann, D. (joint 1st authors), Iglesias, J.E., van der Kouwe, A.J.W.,"
echo "Boyd, E., Reuter, M., Stevens, A., Van Leemput, K., McKee, A., Frosch, M.P., Fischl, B.,"
echo "and Augustinack, J.C., High-resolution magnetic resonance imaging reveals nuclei of the" 
echo "human amygdala: manual segmentation to automatic atlas. Neuroimage 155, 2017, 370-382."
echo "http://doi.org/10.1016/j.neuroimage.2017.04.046"
echo " "

exit 0



