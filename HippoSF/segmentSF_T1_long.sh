#!/bin/tcsh

if (  $#argv < 6 ) then
  echo "Software requires at least 6 arguments"
  echo "segmentSF_T1_long.sh matlabRuntimeDirectory FShomeDirectory subjectDir side baseSubjectName tp1 ... tpN"
  exit 1
endif

# Absolute name of script
set rootdir = `dirname $0`
set SCRIPTPATH = `cd $rootdir && pwd`

# Parameters
set RUNTIME=$1; shift
set FREESURFER_HOME=$1; shift
set SUBJECTDIR = `cd $1 && pwd`; shift
set SIDE=$1; shift
set BASESUBJECTNAME=$1; shift
set RESOLUTION="0.333333333333333333333333333333333333";
set ATLASMESH="$FREESURFER_HOME/average/HippoSF/atlas/AtlasMesh.gz";
set ATLASDUMP="$FREESURFER_HOME/average/HippoSF/atlas/AtlasDump.mgz";
set LUT="$FREESURFER_HOME/average/HippoSF/atlas/compressionLookupTable.txt";
set K1="0.05";
set K2="0.05";
set OPTIMIZER="ConjGrad";
set MRFCONSTANT="0";
set SUFFIX="long.v10";

# command
set cmd="$SCRIPTPATH/run_SegmentSubfieldsT1Longitudinal.sh $RUNTIME $SUBJECTDIR $RESOLUTION $ATLASMESH $ATLASDUMP $LUT $K1 $K2 $SIDE $OPTIMIZER $SUFFIX  ${FREESURFER_HOME}/bin/  $MRFCONSTANT $BASESUBJECTNAME"

foreach file ( $* )
  set cmd="$cmd $file"
end

# echo $cmd

eval $cmd

exit


