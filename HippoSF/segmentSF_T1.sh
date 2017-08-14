#!/bin/tcsh

if (  $#argv < 1 ) then
  echo "Software requires 1 argument:"
  echo ""
  echo "Usage:"
  echo "  segmentSF_T1.sh SubjectName"
  exit 1
endif

checkMCR
if($status) exit 1;

# Parameters
set RUNTIME=$FREESURFER_HOME/MCRv80;
set SUBJECTNAME=$1;
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="ConjGrad";
set MRFCONSTANT="0";
set SUFFIX="v20";

if (! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/wmparc.mgz || \
    ! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/norm.mgz ) then
  echo "ERROR: cannot find norm.mgz or wmparc.mgz for the subject."
  echo "ERROR: Make sure recon-all was run on this subject, to completion."
  exit 1;
endif 

set HSFLOG = ($SUBJECTS_DIR/$SUBJECTNAME/scripts/hippocampal-subfields-T1.log)
rm -f $HSFLOG

# command
foreach SIDE ( left right )
  set cmd="run_segmentSubjectT1_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $SIDE $OPTIMIZER $SUFFIX ${FREESURFER_HOME}/bin/ $MRFCONSTANT"
  eval $cmd |& tee -a $HSFLOG
  if($status) exit 1;
end
 
exit
