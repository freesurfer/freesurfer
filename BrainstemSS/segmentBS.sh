#!/bin/tcsh

if (  $#argv < 1 ) then
  echo "Software requires 1 argument:"
  echo ""
  echo "Usage:"
  echo "  segmentBS.sh SubjectName"
  exit 1
endif

checkMCR
if($status) exit 1;

# Parameters
set RUNTIME=$FREESURFER_HOME/MCRv80;
set SUBJECTNAME=$1;
set RESOLUTION="0.5";
set ATLASMESH="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/BrainstemSS/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="ConjGrad";
set SUFFIX="v10";

# brainstem structures: make sure required files from from default freesurfer processing are in there
if (! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/aseg.mgz || \
    ! -e $SUBJECTS_DIR/$SUBJECTNAME/mri/norm.mgz ) then
  echo "ERROR: cannot find norm.mgz or aseg.mgz for the subject."
  echo "ERROR: Make sure recon-all was run on this subject, to completion."
  exit 1;
endif #checking of aseg.mgz and norm.mgz 

# command
set cmd="run_SegmentSubject.sh $RUNTIME $SUBJECTNAME $SUBJECTS_DIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $OPTIMIZER $SUFFIX ${FREESURFER_HOME}/bin/"

set BSSLOG = ($SUBJECTS_DIR/$SUBJECTNAME/scripts/brainstem-structures.log)
rm -f $BSSLOG
eval $cmd |& tee -a $BSSLOG

exit $?
