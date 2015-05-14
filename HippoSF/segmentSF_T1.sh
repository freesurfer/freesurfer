#!/bin/tcsh

if (  $#argv < 5 ) then
  echo "Software requires 5 arguments"
  echo "segmentSF_T1.sh matlabRuntimeDirectory FShomeDirectory subjectName subjectDir side(must be left or right, in lower case) "
  exit 1
endif

# Absolute name of script
set rootdir = `dirname $0`
set SCRIPTPATH = `cd $rootdir && pwd`

# Parameters
set RUNTIME=$1;
set FREESURFER_HOME=$2;
set SUBJECTNAME=$3;
set SUBJECTDIR = `cd $4 && pwd`
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K="0.05";
set SIDE=$5;
set OPTIMIZER="ConjGrad";
set MRFCONSTANT="0";
set SUFFIX="v1.0";

# command
set cmd="$SCRIPTPATH/run_segmentSubjectT1_autoEstimateAlveusML.sh $RUNTIME $SUBJECTNAME $SUBJECTDIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $SIDE $OPTIMIZER $SUFFIX ${FREESURFER_HOME}/bin/ $MRFCONSTANT"

# echo $cmd

eval $cmd

exit
