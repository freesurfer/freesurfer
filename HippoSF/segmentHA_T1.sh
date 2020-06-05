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
if( $#argv == 0 || $#argv > 2) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   segmentHA_T1.sh SUBJECT_ID [SUBJECT_DIR]" 
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentHA_T1.sh --help"  
  echo " "
  exit 1
endif

# If requesting help
if( $1 == "--help") then
  echo " "
  echo "SEGMENTATION OF HIPPOCAMPAL SUBFIELDS AND NUCLEI OF THE AMYGDALA (T1)"
  echo " "
  echo "Use this script to segment the hippocampal subfields and nuclei of the amygdala"
  echo "from the main T1 scan used in the recon-all stream. It supports standard (1mm)"
  echo "and high resolution T1 scans."
  echo " "
  echo "To use this module, you first need to process the input scan with the main "
  echo "FreeSurfer pipeline (recon-all). If you want to analyze high-resolution scans"
  echo "(i.e., voxel size smaller than 1mm) at their native resolution, they must be"
  echo "processed with the recon-all flag -cm. Then, you can run the following command"
  echo "to obtain the segmentation: "
  echo " "
  echo "   segmentHA_T1.sh SUBJECT_ID [SUBJECT_DIR]"
  echo " "
  echo "   (the argument [SUBJECT_DIR] is only necessary if the"
  echo "    environment variable SUBJECTS_DIR has not been set"
  echo "    or if you want to override it)"
  echo " "
  echo "The segmentation method is described in [1], and the atlases are described in [1]"
  echo "(hippocampus) and [2] (amygdala):"
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "https://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalSubfieldsAndNucleiOfAmygdala"
  echo " "
  echo "[1] Iglesias, J.E., Augustinack, J.C., Nguyen, K., Player, C.M., Player, A., Wright,"
  echo "M., Roy, N., Frosch, M.P., McKee, A.C., Wald, L.L., Fischl, B., and Van Leemput, K.,"
  echo "A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution"
  echo "MRI: Application to adaptive segmentation of in vivo MRI.  Neuroimage 115, 2015, 117-137." 
  echo "http://dx.doi.org/10.1016/j.neuroimage.2015.04.042"
  echo " "
  echo "[2] Saygin, Z.M. & Kliemann, D. (joint 1st authors), Iglesias, J.E., van der Kouwe, A.J.W.,"
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
  set SUBJECTS_DIR = $2
endif

# Set name of subject
set SUBJECTNAME = $1

# Error if subject directory does not exist
if (! -d $SUBJECTS_DIR ) then
  echo " "
  echo "Subjects directory:"
  echo "   $SUBJECTS_DIR"
  echo "does not exist"
  echo " "
  exit 1
endif

# Error if directory of specific subject does not exist
if (! -d $SUBJECTS_DIR/$SUBJECTNAME ) then
  echo " "
  echo "Directory of subject to process:"
  echo "   $SUBJECTS_DIR/$SUBJECTNAME"
  echo "does not exist"
  echo " "
  exit 1
endif

# Error if subject not processed (far enough)
if (! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/wmparc.mgz || \
    ! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/norm.mgz || \
    ! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/transforms/talairach.xfm ) then
  echo " "
  echo "Cannot find wmparc.mgz or norm.mgz or talairach.xfm for the subject."
  echo "Has the subject been procesed with recon-all?"
  echo " "
  exit 1;
endif 

# Make sure that the (T1) hippocampal subfields are not running already for this subject
set IsRunningFile = ${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/IsRunningHPsubT1.lh+rh
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that T1 Hippocampal Subfields is already running"
  echo "for this subject based on the presence of IsRunningHPsub.lh+rh"
  echo "It could also be that T1 Hippocampal Subfields was running at one"
  echo "point but died in an unexpected way. If it is the case that there"
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

# If not explicitly specfied, set to 1
if($?ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS == 0) then
  setenv ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS 1
endif

# If everything is in place, let's do it! First, we create the IsRunning file
echo "------------------------------" > $IsRunningFile
echo "SUBJECT $SUBJECTNAME" >> $IsRunningFile 
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
set RUNTIME="$FREESURFER_HOME/MCRv84/";
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="L-BFGS";
set MRFCONSTANT="0";
set SUFFIX="v21";

# Now the real job
set hippohemilist=(left right)
set HSFLOG = (${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/hippocampal-subfields-T1.log)
rm -f $HSFLOG

echo "------------------------------" > $HSFLOG
echo "USER $user"      >> $HSFLOG
echo "HOST `hostname`" >> $HSFLOG
echo "PROCESSID $$ "   >> $HSFLOG
echo "PROCESSOR `uname -m`" >> $HSFLOG
echo "OS `uname -s`"       >> $HSFLOG
echo "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS " >> $HSFLOG
uname -a         >> $HSFLOG
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $HSFLOG
endif
echo "------------------------------" >> $HSFLOG
echo " " >> $HSFLOG
cat $FREESURFER_HOME/build-stamp.txt  >> $HSFLOG
echo " " >> $HSFLOG
echo "setenv SUBJECTS_DIR $SUBJECTS_DIR"  >> $HSFLOG
echo "cd `pwd`"   >> $HSFLOG
echo $0 $argv  >> $HSFLOG
echo ""  >> $HSFLOG

foreach hemi ($hippohemilist)
  echo "#--------------------------------------------" \
    |& tee -a $HSFLOG
  echo "#@# Hippocampal Subfields processing (T1) $hemi `date`" \
    |& tee -a $HSFLOG

  # command
  set cmd="run_segmentSubjectT1_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $hemi $OPTIMIZER $SUFFIX '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/'  $MRFCONSTANT"

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
    echo "T1 hippocampal subfields exited with ERRORS at `date`" \
    |& tee -a $HSFLOG
    echo "" |& tee -a $HSFLOG
    echo "For more details, see the log file $HSFLOG"
    echo ""
    rm -f $IsRunningFile
    exit 1;
  endif
end

# Convert the txt files into a stats file so that asegstats2table can
# be run Note: the number of voxels is set to 0 and there is no info
# about intensity. The only useful info is the volume in mm and the
# structure name. The segmentation IDs also do not mean anything. 
# Could run mri_segstats instead, but the volumes would not include
# partial volume correction.
foreach hemi (lh rh)
  set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/$hemi.hippoSfVolumes-T1.$SUFFIX.txt
  set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/hipposubfields.$hemi.T1.$SUFFIX.stats
  echo "# Hippocampal subfield volumes as created by segmentHA_T1.sh" > $stats
  cat $txt | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats
  set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/$hemi.amygNucVolumes-T1.$SUFFIX.txt
  set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/amygdalar-nuclei.$hemi.T1.$SUFFIX.stats
  echo "# Amygdala nuclei volumes as created by segmentHA_T1.sh" > $stats
  cat $txt | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats
end
 
# All done!
rm -f $IsRunningFile

echo " "
echo "All done!"
echo " "
echo "If you have used results from this software for a publication, please cite:"
echo " "
echo "Iglesias, J.E., Augustinack, J.C., Nguyen, K., Player, C.M., Player, A., Wright,"
echo "M., Roy, N., Frosch, M.P., McKee, A.C., Wald, L.L., Fischl, B., and Van Leemput, K.,"
echo "A computational atlas of the hippocampal formation using ex vivo, ultra-high resolution"
echo "MRI: Application to adaptive segmentation of in vivo MRI.  Neuroimage 115, 2015, 117-137." 
echo "http://dx.doi.org/10.1016/j.neuroimage.2015.04.042"
echo " "
echo "In addition, if you have used the segmentation of the nuclei of the amygdala, please also cite:"
echo ""
echo "Saygin, Z.M. & Kliemann, D. (joint 1st authors), Iglesias, J.E., van der Kouwe, A.J.W.,"
echo "Boyd, E., Reuter, M., Stevens, A., Van Leemput, K., McKee, A., Frosch, M.P., Fischl, B.,"
echo "and Augustinack, J.C., High-resolution magnetic resonance imaging reveals nuclei of the" 
echo "human amygdala: manual segmentation to automatic atlas. Neuroimage 155, 2017, 370-382."
echo "http://doi.org/10.1016/j.neuroimage.2017.04.046"
echo " "

exit 0



