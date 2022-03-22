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
  echo "SEGMENTATION OF THALAMIC NUCLEI"
  echo " "
  echo "Use this script to segment the thalamic nuclei from either the main"
  echo "T1 scan used in the recon-all stream, or an additional volume with"
  echo "a different type of MRI contrast (e.g., FGATIR, which yields great"
  echo "contrast in the thalamus)."
  echo " "
  echo "To use this module, you first need to process the input scan with the main "
  echo "FreeSurfer pipeline (recon-all). If you want to analyze high-resolution scans"
  echo "(i.e., voxel size smaller than 1mm) at their native resolution, they must be"
  echo "processed with the recon-all flag -cm. Then, you can run the following command"
  echo "to obtain the thalamic segmentation: "
  echo " "
  echo "segmentThalamus.sh SUBJECT_ID [SUBJECT_DIR] [ADDITIONAL_VOL   ANALYSIS_ID   BBREGISTER_MODE] "
  echo "  "
  echo "  SUBJECT_ID: FreeSurfer subject name, e.g., bert "
  echo '  SUBJECT_DIR: FreeSurfer subjects directory, typically \$SUBJECTS_DIR '
  echo "               (the argument [SUBJECT_DIR] is only necessary if the"
  echo "                environment variable SUBJECTS_DIR has not been set"
  echo "                or if you want to override it)"
  echo "  "
  echo "  Optional arguments for segmentation of other MRI contrast (additional volume) "
  echo "  "
  echo "  ADDITIONAL_VOL: scan to segment (supported formats: .nii, .nii.gz, .mgz, .mgh) "
  echo '  SUBJECT_DIR: FreeSurfer subjects directory, typically $SUBJECTS_DIR'
  echo "  ANALYSIS_ID: a user-defined string describing the type of sequence use. For example: T2, FGATIR, FLAIR ..."
  echo "  BBREGISTER_MODE: contrast of additional scan, used by bbregister in the registration. Must be 't1', 't2', or 'none': "
  echo "                  t1: Assume t1-like contrast, ie, White Matter brighter than Grey Matter "
  echo "                  t2: Assume t2-like contrast, ie, Gray Matter brighter than White Matter "
  echo "                  none: Assume that input scan and FreeSurfer's main T1 are already registered "
  echo "                        (skips registration, supports any MRI contrast) "
  echo " "
  echo "See further information, including how to visualize the results, at:"
  echo " "
  echo "http://surfer.nmr.mgh.harvard.edu/fswiki/ThalamicNuclei"
  echo " "
  echo "The segmentation method is described in the following publication:"
  echo " "
  echo "A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology "
  echo "Iglesias, J.E., Insausti, R., Lerma-Usabiaga, G., Bocchetta, M.,"
  echo "Van Leemput, K., Greve, D., van der Kouwe, A., Caballero-Gaudes, C., "
  echo "Paz-Alonso, P. NeuroImage (in press)."
  echo "Preprint available at arXiv.org:  https://arxiv.org/abs/1806.08634" 
  echo " "
  exit 0
endif

# If wrong number of arguments given
if( $#argv != 1  && $#argv != 2 && $#argv != 5 ) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   segmentThalamicNuclei.sh SUBJECT_ID [SUBJECT_DIR] [ADDITIONAL_VOL   ANALYSIS_ID   BBREGISTER_MODE] " 
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentThalamicNuclei.sh --help"  
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

# Set additional parameters, if necessary
set ANALYSISID="mainFreeSurferT1";
if ($#argv > 2) then
  set ADDVOL="`getfullpath  $3`";
  set ANALYSISID="$4";
  set BBREGMODE="$5";
  set DOBIASFIELDCORR="1";
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

# Make sure that the thalamic nuclei are not running already for this subject
set IsRunningFile = ${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/IsRunningThalamicNuclei_${ANALYSISID}
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that thalamic nuclei segmentation is already running"
  echo "for this subject based on the presence of IsRunningThalamicNuclei_${ANALYSISID}"
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
set RESOLUTION="0.5";
set ATLASMESH="$FREESURFER_HOME/average/ThalamicNuclei/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/ThalamicNuclei/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/ThalamicNuclei/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="L-BFGS";
set SUFFIX="v13";
set USETWOCOMPS="1";
set MRFCONSTANT="0";


# Now the real job
set THNUCLOG = (${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/thalamic-nuclei-${ANALYSISID}.log)
rm -f $THNUCLOG

echo "------------------------------" > $THNUCLOG
echo "USER $user"      >> $THNUCLOG
echo "HOST `hostname`" >> $THNUCLOG
echo "PROCESSID $$ "   >> $THNUCLOG
echo "PROCESSOR `uname -m`" >> $THNUCLOG
echo "OS `uname -s`"       >> $THNUCLOG
echo "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS " >> $THNUCLOG
uname -a         >> $THNUCLOG
if($?PBS_JOBID) then
  echo "pbsjob $PBS_JOBID"  >> $THNUCLOG
endif
echo "------------------------------" >> $THNUCLOG
cat $FREESURFER_HOME/build-stamp.txt  >> $THNUCLOG
echo " " >> $THNUCLOG
echo "setenv SUBJECTS_DIR $SUBJECTS_DIR"  >> $THNUCLOG
echo "cd `pwd`"   >> $THNUCLOG
echo $0 $argv  >> $THNUCLOG
echo ""  >> $THNUCLOG

echo "#--------------------------------------------" \
  |& tee -a $THNUCLOG
echo "#@# Thalamic Nuclei processing `date`" \
  |& tee -a $THNUCLOG

# command
set cmd="run_SegmentThalamicNuclei.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $OPTIMIZER $SUFFIX '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/' $USETWOCOMPS  $MRFCONSTANT"
if ($#argv > 2) then
  set cmd="$cmd $ADDVOL $ANALYSISID $DOBIASFIELDCORR $BBREGMODE";
endif

fs_time ls >& /dev/null
if ($status) then
  $cmd |& tee -a $THNUCLOG 
  set returnVal=$status
else
  fs_time $cmd |& tee -a $THNUCLOG
  set returnVal=$status
endif
echo $returnVal


if ($returnVal) then
  uname -a | tee -a $THNUCLOG
  echo "" |& tee -a $THNUCLOG
  echo "Thalamic nuclei exited with ERRORS at `date`" \
  |& tee -a $THNUCLOG
  echo "" |& tee -a $THNUCLOG
  echo "For more details, see the log file $THNUCLOG"
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
if("$ANALYSISID" == "mainFreeSurferT1") then
  set suffix2 = $SUFFIX.T1
else
  set suffix2 = $SUFFIX.$ANALYSISID
endif
set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/ThalamicNuclei.$suffix2.volumes.txt
# Divide the stats into left and right. The sorting is done because
# the Left and Right nuclei are not ordered the same in the txt file
set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/thalamic-nuclei.lh.$suffix2.stats
echo "# Left Thalamic nuclei volume statistics as created by segmentThalamicNuclei.sh" > $stats
grep Left $txt | sed 's/Left-//g' | sort -k 1 | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats
set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/thalamic-nuclei.rh.$suffix2.stats
echo "# Right Thalamic nuclei volume statistics as created by segmentThalamicNuclei.sh" > $stats
grep Right $txt | sed 's/Right-//g' | sort -k 1 | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats

# All done!
rm -f $IsRunningFile

echo " "
echo "All done!"
echo " "
echo "If you have used results from this software for a publication, please cite:"
echo " "
echo "A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology "
echo "Iglesias, J.E., Insausti, R., Lerma-Usabiaga, G., Bocchetta, M.,"
echo "Van Leemput, K., Greve, D., van der Kouwe, A., Caballero-Gaudes, C., "
echo "Paz-Alonso, P. NeuroImage (in press)."
echo "Preprint available at arXiv.org:  https://arxiv.org/abs/1806.08634"
echo " "

exit 0

