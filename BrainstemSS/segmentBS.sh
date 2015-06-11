#!/bin/tcsh



if (  $#argv < 4 ) then
  echo "Software requires 4 arguments"
  echo "segmentSFnewT1.sh matlabRuntimeDirectory FShomeDirectory subjectName subjectDir  "
  exit 1
endif


# Absolute name of script
set rootdir = `dirname $0`
set SCRIPTPATH = `cd $rootdir && pwd`

## Environment variables and FreeSurfer
#source ${FREESURFER_HOME}/SetUpFreeSurfer.sh
#PATH=$PATH:$SCRIPTPATH/
#export PATH


# Parameters
set RUNTIME=$1;
set FREESURFER_HOME=$2;
set SUBJECTNAME=$3;
set SUBJECTDIR = `cd $4 && pwd`
set RESOLUTION="0.5";
set ATLASMESH="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/BrainstemSS/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/BrainstemSS/atlas/compressionLookupTable.txt";
set K="0.05";
set OPTIMIZER="ConjGrad";
set SUFFIX="v10";

# command
set cmd="$SCRIPTPATH/run_SegmentSubject.sh $RUNTIME $SUBJECTNAME $SUBJECTDIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K $OPTIMIZER $SUFFIX ${FREESURFER_HOME}/bin/"

# echo $cmd

eval $cmd

exit



