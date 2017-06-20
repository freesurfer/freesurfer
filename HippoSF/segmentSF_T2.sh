#!/bin/tcsh

if (  $#argv < 3 ) then
  echo "Software requires 3 arguments:"
  echo ""
  echo "Usage:"
  echo "  segmentSF_T1T2.sh SubjectName FileNameAdditionalFile AnalysisID"
  exit 1
endif

checkMCR
if($status) exit 1;

# Parameters
set RUNTIME=$FREESURFER_HOME/MCRv80;
set SUBJECTNAME=$1;
set T2VOL=$2
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="ConjGrad";
set MRFCONSTANT="0";
set SUFFIX="v20";
set USERSUFFIX=$3;
set BYPASSBF="1";
set USEWHOLEBRAININHP="0";

if (! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/wmparc.mgz || \
    ! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/norm.mgz ) then
  echo "ERROR: cannot find norm.mgz or wmparc.mgz for the subject."
  echo "ERROR: Make sure recon-all was run on this subject, to completion."
  exit 1;
endif 

set HSFLOG = ($SUBJECTS_DIR/$SUBJECTNAME/scripts/hippocampal-subfields-T2.log)
rm -f $HSFLOG

# command
foreach SIDE ( left right )
# command
  set cmd="$SCRIPTPATH/run_segmentSubjectT2_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $T2VOL $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $SIDE $OPTIMIZER $SUFFIX $USERSUFFIX  ${FREESURFER_HOME}/bin/ $MRFCONSTANT $BYPASSBF $USEWHOLEBRAININHP"
  eval $cmd |& tee -a $HSFLOG
  if($status) exit 1;
end

exit
