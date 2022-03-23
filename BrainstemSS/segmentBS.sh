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
  echo "   segmentBS.sh SUBJECT_ID [SUBJECT_DIR]" 
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentBS.sh --help"  
  echo " "
  exit 1
endif

# If requesting help
if( $1 == "--help") then
  echo " "
  echo "SEGMENTATION OF BRAINSTEM SUBSTRUCTURES"
  echo " "
  echo "Use this script to segment the brainstem substructure (medulla, pons, midbrain,"
  echo "and SCP) from the main T1 scan used in the recon-all stream."
  echo " "
  echo "To use this module, you first need to process the input scan with the main "
  echo "FreeSurfer pipeline (recon-all). If you want to analyze high-resolution scans"
  echo "(i.e., voxel size smaller than 1mm) at their native resolution, they must be"
  echo "processed with the recon-all flag -cm. Then, you can run the following command"
  echo "to obtain the segmentation: "
  echo " "
  echo "   segmentBS.sh SUBJECT_ID [SUBJECT_DIR]"
  echo " "
  echo "   (the argument [SUBJECT_DIR] is only necessary if the"
  echo "    environment variable SUBJECTS_DIR has not been set"
  echo "    or if you want to override it)"
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "http://surfer.nmr.mgh.harvard.edu/fswiki/BrainstemSubstructures"
  echo " "
  echo "The segmentation method is described in the following publication:"
  echo " "
  echo "Iglesias, J.E., Van Leemput, K., Bhatt, P., Casillas, C., Dutt, S., Schuff, N.,"
  echo "Truran-Sacrey, D., Boxer, A., and Fischl, B., Bayesian segmentation of brainstem "
  echo "structures in MRI. Neuroimage 113, 2015, 184-195."
  echo "http://dx.doi.org/10.1016/j.neuroimage.2015.02.065"
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
set IsRunningFile = ${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/IsRunningBSsubst
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that Brainstem Substructures is already running"
  echo "for this subject based on the presence of IsRunningBSsubst"
  echo "It could also be that Brainstem Substructures was running at one"
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
set RUNTIME="$FREESURFER_HOME/MCRv97/";
set RESOLUTION="0.5";
set ATLASMESH="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/BrainstemSS/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="L-BFGS";
set MRFCONSTANT="0";
set SUFFIX="v13";


# Now the real job
set BSSLOG = (${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/brainstem-substructures-T1.log)
rm -f $BSSLOG

echo "------------------------------" > $BSSLOG
echo "USER $user"      >> $BSSLOG
echo "HOST `hostname`" >> $BSSLOG
echo "PROCESSID $$ "   >> $BSSLOG
echo "PROCESSOR `uname -m`" >> $BSSLOG
echo "OS `uname -s`"       >> $BSSLOG
echo "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS " >> $BSSLOG
uname -a         >> $BSSLOG
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $BSSLOG
endif
echo "------------------------------" >> $BSSLOG
echo " " >> $BSSLOG
cat $FREESURFER_HOME/build-stamp.txt  >> $BSSLOG
echo " " >> $BSSLOG
echo "setenv SUBJECTS_DIR $SUBJECTS_DIR"  >> $BSSLOG
echo "cd `pwd`"   >> $BSSLOG
echo $0 $argv  >> $BSSLOG
echo ""  >> $BSSLOG

echo "#--------------------------------------------" \
  |& tee -a $BSSLOG
echo "#@# Brainstem Substructures processing `date`" \
  |& tee -a $BSSLOG

# command
set cmd="run_SegmentSubject.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $OPTIMIZER $SUFFIX '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/'"

fs_time ls >& /dev/null
if ($status) then
  $cmd |& tee -a $BSSLOG 
  set returnVal=$status
else
  fs_time $cmd |& tee -a $BSSLOG
  set returnVal=$status
endif

if ($returnVal) then
  uname -a | tee -a $BSSLOG
  echo "" |& tee -a $BSSLOG
  echo "T1 Brainstem Substructures exited with ERRORS at `date`" \
  |& tee -a $BSSLOG
  echo "" |& tee -a $BSSLOG
  echo "For more details, see the log file $BSSLOG"
  echo ""
  rm -f $IsRunningFile
  exit 1;
endif

# Convert the txt file into a stats file so that asegstats2table can
# be run Note: the number of voxels is set to 0 and there is no info
# about intensity. The only useful info is the volume in mm and the
# structure name. The segmentation IDs also do not mean anything. 
# Could run mri_segstats instead, but the volumes would not include
# partial volume correction.
set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/brainstemSsVolumes.$SUFFIX.txt
set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/brainstem.$SUFFIX.stats
echo "# Brainstem structure volumes as created by segmentBS.sh" > $stats
cat $txt | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats

# All done!
rm -f $IsRunningFile

echo " "
echo "All done!"
echo " "
echo "If you have used results from this software for a publication, please cite:"
echo " "
  echo "Iglesias, J.E., Van Leemput, K., Bhatt, P., Casillas, C., Dutt, S., Schuff, N.,"
  echo "Truran-Sacrey, D., Boxer, A., and Fischl, B., Bayesian segmentation of brainstem "
  echo "structures in MRI. Neuroimage 113, 2015, 184-195."
  echo "http://dx.doi.org/10.1016/j.neuroimage.2015.02.065"
echo " "

exit 0
