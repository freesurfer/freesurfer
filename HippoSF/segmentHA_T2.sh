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
if( $#argv == 0 || $#argv == 2  || $#argv == 3 || $#argv > 5) then
  echo " "
  echo "Usage: "
  echo " "
  echo "   segmentHA_T2.sh  SUBJECT_ID  FILE_ADDITIONAL_SCAN   ANALYSIS_ID  USE_T1  [SUBJECT_DIR]" 
  echo " "
  echo "Or, for help"
  echo " "
  echo "   segmentHA_T2.sh --help"  
  echo " "
  exit 1
endif

# If requesting help
if( $1 == "--help") then
  echo " "
  echo "SEGMENTATION OF HIPPOCAMPAL SUBFIELDS AND NUCLEI OF THE AMYGDALA (T2 or T1+T2)"
  echo " "
  echo "Use this script for segmentation of the hippocampal subfields and nuclei of"
  echo "the amygdala with an additional scan, which is typically (but not"
  echo "necessarily) a higher resolution T2 scan covering the hippocampal formation."
  echo "This script has two modes of operation: using the intensities of the "
  echo "additional scan alone, or combining them with the intensities of the  main "
  echo "T1 input scan from recon-all (multispectral segmentation). The former is "
  echo "recommended when the additional scan is of higher resolution than the T1 in"
  echo "every direction and covers the whole hippocampal ROI; if this is not the case,"
  echo "using both scans simultaneously normally produces better results."
  echo " "
  echo "To use this module, you first need to process the main T1 input scan with the "
  echo "main FreeSurfer pipeline (recon-all). Bear in mind that, if you want to "
  echo "analyze high-resolution scans (i.e., voxel size smaller than 1mm) at their"
  echo "native resolution, they must be processed with the recon-all flag -cm. Then,"
  echo "you can run the following command to obtain the segmentation: "
  echo " "
  echo "   segmentHA_T2.sh  SUBJECT_ID  FILE_ADDITIONAL_SCAN   ANALYSIS_ID  USE_T1  [SUBJECT_DIR]"
  echo " "
  echo "   (the argument [SUBJECT_DIR] is only necessary if the"
  echo "    environment variable has not been set)"
  echo " "
  echo "FILE_ADDITIONAL_SCAN is the additional scan in Nifti (.nii/.nii.gz) or FreeSurfer"
  echo "format (.mgh/.mgz)."
  echo " "
  echo "ANALYSIS_ID is a user defined identifier that makes it possible to run different"
  echo "analysis with different types of additional scans, e.g., a T2 and proton density"
  echo "weighted volume (see further information on the website below). Please note that it". 
  echo "is not possible to two instances of run segmentHA_T2.sh in parallel with USE_T1=0 "
  echo "and USE_T1=1  (segmentHA_T1.sh or different analysis IDs is fine)."
  echo " "
  echo "USE_T1: set to 1 to use the main T1 from recon-all in the segmentation (multispectral "
  echo "mode). If set to 0, only the intensities of the additional scan are used".
  echo " "
  echo "SUBJECT_DIR (optional argument): FreeSurfer subject directory (overrides environment"
  echo "variable SUBJECTS_DIR)."
  echo " "
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
if ($#argv == 4) then
  if (! $?SUBJECTS_DIR)  then
    echo " "
    echo "SUBJECTS_DIR variable does not exist"
    echo "Please define it or provide subjects directory as second input"
    echo " "
    exit 1
  endif
endif

# Error if SUBJECTS_DIR (the environemnt variable) is empty 
if ($#argv == 4) then
  if ( $SUBJECTS_DIR == "" ) then
    echo " "
    echo "SUBJECTS_DIR variable is empty"
    echo "Please redefine it or provide subjects directory as second input"
    echo " "
    exit 1
  endif
endif

# If SUBJECTS_DIR is provided, just set it
if ($#argv == 5) then
  set SUBJECTS_DIR = `getfullpath  $5`
endif

# Set name of subject
set SUBJECTNAME = $1
set T2VOL = `getfullpath  $2`
set ANALYSISID = $3
set USET1 = $4

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

# Error if additional volume does not exist
if(! -e "$T2VOL") then
  echo "ERROR: cannot find additional scan $T2VOL"
  echo " "
  exit 1;
endif

# Error if additional volume is a directory
if(-d "$T2VOL") then
  echo "ERROR: additional scan $T2VOL cannot be a directory"
  echo " "
  exit 1;
endif

# Error if additional volume is not readable
if(! -r "$T2VOL") then
  echo "ERROR: additional scan $T2VOL is not readable"
  echo " "
  exit 1;
endif

# Error if USET1 is not 0 or 1
if (! ($USET1 =~ {"0","1"}) ) then
  echo "ERROR: USE_T1 must be 0 or 1"
  echo " "
  exit 1;
endif

# Make sure that the T2 (+T1) hippocampal subfields are not running already for this subject
# with the same analysis ID (other analysis ID is fine)
set IsRunningFile = ${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/IsRunningHPsubT2.${ANALYSISID}.lh+rh
if(-e $IsRunningFile) then
  echo ""
  echo "It appears that T2 (+T1) Hippocampal Subfields is already running"
  echo "for this subject with the same analysis ID, based on the presence of"
  echo "IsRunningHPsubT2.${ANALYSISID}.lh+rh"
  echo "It could also be that T2 Hippocampal Subfields was running"
  echo "at one point but died in an unexpected way. If it is the case that "
  echo "there is a process running, you can kill it and start over or just"
  echo "let it run. If the process has died, you should type:"
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
echo "SUBJECT $SUBJECTNAME" >> $IsRunningFile 
echo "ANALYSIS_ID $ANALYSISID" >> $IsRunningFile 
echo "USE_T1 $USET1" >> $IsRunningFile 
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
set KT2="0.05";
set KT1T2="0.10";
set OPTIMIZER="L-BFGS";
set MRFCONSTANT="0";
set SUFFIX="v21";
set BYPASSBF="1";
set USEWHOLEBRAININHP="0";

# Now the real job
set hippohemilist=(left right)
set HSFLOG = (${SUBJECTS_DIR}/${SUBJECTNAME}/scripts/hippocampal-subfields-T2.${ANALYSISID}.log)
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
  

  # command
  if ( $USET1 == "1") then
    echo "#@# Hippocampal Subfields processing (T1+T2) $hemi `date`"  |& tee -a $HSFLOG
    set cmd="run_segmentSubjectT1T2_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $T2VOL $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $KT1T2 $hemi $OPTIMIZER $SUFFIX $ANALYSISID '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/' $MRFCONSTANT $BYPASSBF $USEWHOLEBRAININHP"
  else
    echo "#@# Hippocampal Subfields processing (T2) $hemi `date`"  |& tee -a $HSFLOG
    set cmd="run_segmentSubjectT2_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $T2VOL $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $KT2 $hemi $OPTIMIZER $SUFFIX $ANALYSISID  '${FREESURFER_HOME}/bin/fs_run_from_mcr ${FREESURFER_HOME}/bin/' $MRFCONSTANT $BYPASSBF $USEWHOLEBRAININHP"
  endif

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
    echo "T2 hippocampal subfields exited with ERRORS at `date`" \
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


  set suffix2 = $SUFFIX.
foreach hemi (lh rh)
  set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/$hemi.hippoSfVolumes-T1-$ANALYSISID.$SUFFIX.txt
  set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/hipposubfields.$hemi.T2.$SUFFIX.$ANALYSISID.stats
  echo "# Hippocampal subfield volumes as created by segmentHA_T2.sh" > $stats
  cat $txt | awk '{print NR" "NR"  0 "$2" "$1}' >> $stats
  set txt=$SUBJECTS_DIR/$SUBJECTNAME/mri/$hemi.amygNucVolumes-T1-$ANALYSISID.$SUFFIX.txt
  set stats=$SUBJECTS_DIR/$SUBJECTNAME/stats/amygdalar-nuclei.$hemi.T2.$SUFFIX.$ANALYSISID.stats
  echo "# Amygdala nuclei volumes as created by segmentHA_T2.sh" > $stats
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






