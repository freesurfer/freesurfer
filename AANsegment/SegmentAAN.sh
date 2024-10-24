#! /bin/tcsh -f

set tcsh61706 = (`tcsh --version | grep "6\.17\.06"`)
if ("$tcsh61706" != "") then
  echo ""
  echo "WARNING: tcsh v6.17.06 has an exit code bug! Please update tcsh!"
  echo ""
  # workaround to force expected behavior:
  set anyerror
endif

# very first thing: check MCR installation
checkMCR.sh
if ($status) then
    exit 1
endif



# If requesting help
if( $1 == "--help") then
  echo " "
  echo "SEGMENTATION OF ASCENDING AROUSAL NETWORK NUCLEI"
  echo " "
  echo "This script segments AAN nuclei from the main"
  echo "T1 scan used in the recon-all stream. Additional contrasts which"
  echo "were tested and seem to work well are T2, low-b (dMRI), MD (dMRI)"
  echo "and SWI."
  echo " "
  echo "To use this module, you first need to process the input scan with the main "
  echo "FreeSurfer pipeline (recon-all). If you want to analyze high-resolution scans"
  echo "(i.e., voxel size smaller than 1mm) at their native resolution, they must be"
  echo "processed with the recon-all flag -cm. Then, you can run the following command"
  echo "You can also side-step the recon-all pipeline (at your own risk!) and create"
  echo "a fake recon-all mri dir, with your image of interest set as the T1. You just"
  echo "need a corresponding ASEG file."
  echo "to obtain the AAN segmentations: "
  echo " "
  echo "segmentAAN.sh SUBJECT_ID [SUBJECT_DIR] "
  echo "  "
  echo "  SUBJECT_ID: FreeSurfer subject name, e.g., bert "
  echo '  SUBJECT_DIR: FreeSurfer subjects directory, typically \$SUBJECTS_DIR '
  echo "               (the argument [SUBJECT_DIR] is only necessary if the"
  echo "                environment variable SUBJECTS_DIR has not been set"
  echo "                or if you want to override it)"
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "http://surfer.nmr.mgh.harvard.edu/fswiki/AANsegment"
  echo " "
  echo "Details on the segmentation method are described in the following preprint:"
  echo " "
  echo "Histology-guided MRI segmentation of brainstem nuclei critical to consciousness"
  echo "Olchanyi, M.D., Augustinack, J., Haynes, R., Lewis, L.D.,"
  echo "Cicero, N., Destrieux, C., Folkerth, R., Kinney, H.C., "
  echo "Fischl, B., Iglesias, J.E., Edlow, B. medRxiv."
  echo "Preprint available at medRxiv.org: https://www.medrxiv.org/content/10.1101/2024.09.26.24314117v1"
  echo " "
  exit 0
endif

# If wrong number of arguments given
if( $#argv != 1  && $#argv != 2 ) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   segmentAAN.sh SUBJECT_ID [SUBJECT_DIR]  "
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentAAN.sh --help"
  echo " "
  exit 1
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

# Warning if scripts directory does not exist
if (! -d $SUBJECTS_DIR/$SUBJECTNAME/scripts ) then
  echo " "
  echo "WARNING: ./scripts Directory in subject to process:"
  echo "   $SUBJECTS_DIR/$SUBJECTNAME"
  echo "does not exist."
  echo "Check if $SUBJECTNAME has been fully processed with recon-all!"
  echo " "
endif

# Error if subject not processed
if (! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/wmparc.mgz || \
    ! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/norm.mgz || \
    ! -e ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/transforms/talairach.xfm ) then
  echo " "
  echo "Cannot find wmparc.mgz or norm.mgz or talairach.xfm for the subject."
  echo "Has the subject been procesed with recon-all?"
  echo " "
  exit 1;
endif

# Make sure that the AAN nuclei are not running already for this subject
set IsRunningFile = ${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/IsRunningAANnuclei
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that AAN nuclei segmentation is already running"
  echo "for this subject based on the presence of IsRunningAANnuclei"
  echo "It could also be that the module was running at one"
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
  echo Setting ITK threads to 1
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
set RUNTIME="$FREESURFER_HOME/MCRv97/";
set RESOLUTION="0.375";
set ATLASMESH="$FREESURFER_HOME/average/AAN/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/AAN/atlas/targetWorkingres.mgz";
set ATLASDUMPWHOLE="$FREESURFER_HOME/average/AAN/atlas/targetReg.mgz";
set LUT="$FREESURFER_HOME/average/AAN/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="L-BFGS";
set SUFFIX="v10";


# Now the real job
set AANNUCLOG = (${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/AAAN-nuclei.log)
rm -f $AANNUCLOG

echo "------------------------------" > $AANNUCLOG
echo "USER $user"      >> $AANNUCLOG
echo "HOST `hostname`" >> $AANNUCLOG
echo "PROCESSID $$ "   >> $AANNUCLOG
echo "PROCESSOR `uname -m`" >> $AANNUCLOG
echo "OS `uname -s`"       >> $AANNUCLOG
echo "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS " >> $AANNUCLOG
uname -a         >> $AANNUCLOG
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $AANNUCLOG
endif
echo "------------------------------" >> $AANNUCLOG
cat $FREESURFER_HOME/build-stamp.txt  >> $AANNUCLOG
echo " " >> $AANNUCLOG
echo "setenv SUBJECTS_DIR $SUBJECTS_DIR"  >> $AANNUCLOG
echo "cd `pwd`"   >> $AANNUCLOG
echo $0 $argv  >> $AANNUCLOG
echo ""  >> $AANNUCLOG

echo "#--------------------------------------------" \
  |& tee -a $AANNUCLOG
echo "#@# AAN Nuclei processing `date`" \
  |& tee -a $AANNUCLOG

# command
set cmd="run_segmentNuclei.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $ATLASDUMPWHOLE $LUT $K $OPTIMIZER $SUFFIX ${FREESURFER_HOME}/bin/"

fs_time ls >& /dev/null
if ($status) then
  $cmd |& tee -a $AANNUCLOG
  set returnVal=$status
else
  fs_time $cmd |& tee -a $AANNUCLOG
  set returnVal=$status
endif
echo $returnVal

# TEMPORARY: convert to .stats file under /stats
mkdir -p ${SUBJECTS_DIR}/${SUBJECTNAME}/stats
mv ${SUBJECTS_DIR}/${SUBJECTNAME}/mri/arousalNetworkVolumes.${SUFFIX}.txt ${SUBJECTS_DIR}/${SUBJECTNAME}/stats/arousalNetworkVolumes.${SUFFIX}.stats

if ($returnVal) then
  uname -a | tee -a $AANNUCLOG
  echo "" |& tee -a $AANNUCLOG
  echo "AAN nuclei exited with ERRORS at `date`" \
  |& tee -a $AANNUCLOG
  echo "" |& tee -a $AANNUCLOG
  echo "For more details, see the log file $AANNUCLOG"
  echo ""
  rm -f $IsRunningFile
  exit 1;
endif

# All done!
rm -f $IsRunningFile



echo " "
echo "All done!"
echo " "
echo "To visualize the outputs, run: "
echo "freeview -v $SUBJECTS_DIR/$SUBJECTNAME/T1.mgz -v $SUBJECTS_DIR/$SUBJECTNAME/arousalNetworkLabels.$SUFFIX.mgz:colormap=lut:lut=$FREESURFER_HOME/average/AAN/atlas/freeview.lut.txt"
echo " "
echo "If you have used results from this software for a publication, please cite:"
echo " "
echo "Histology-guided MRI segmentation of brainstem nuclei critical to consciousness"
echo "Olchanyi, M.D., Augustinack, J., Haynes, R., Lewis, L.D.,"
echo "Cicero, N., Destrieux, C., Folkerth, R., Kinney, H.C., "
echo "Fischl, B., Iglesias, J.E., Edlow, B. medRxiv."
echo "Preprint available at medRxiv.org: https://www.medrxiv.org/content/10.1101/2024.09.26.24314117v1"
echo " "

exit 0
