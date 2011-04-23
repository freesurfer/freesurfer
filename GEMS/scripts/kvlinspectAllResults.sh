#!/bin/bash


# Which parts are you going to do
doLeft=true;
doRight=true;
doAffineRegistration=true;
doDeformableRegistration=true;
doNonPartialVolumedSegmentation=true;
doPartialVolumedSegmentation=false;


# Parse input
if [ $# -lt 1 ]; then
  echo
  echo "Usage: $0 resultsDirectory [ startIndex=1 endIndex=end ]"
  echo
  echo "    e.g.:"
  echo
  echo "       $0 alexResults 1 4"
  echo
  exit -1
fi
resultsDirectory=$1


# Define a function that will show what we're doing, and
# exits if something goes wrong. 
function doIt {

  command=$1

  echo $command
  eval $command

  if [ $? != 0 ]; then
    echo "failed to do $command"
  exit -1
fi

}


doIt "cd $resultsDirectory"

startIndex=1
if [ $# -ge 2 ]; then
  startIndex=$2
fi

subjectNames=(`ls -d */`) # list only directories
if [ $? != 0 ]; then
  echo "failed to list file names in $resultsDirectory"
  exit -1
fi
numberOfSubjects=${#subjectNames[*]}
#echo "numberOfSubjects: $numberOfSubjects"

endIndex=$numberOfSubjects
if [ $# -ge 3 ]; then
  endIndex=$3
  if [ $endIndex -gt $numberOfSubjects ]; then
    echo "endIndex ($endIndex) can't be greater than the number of subjects ($numberOfSubjects)"
    exit -1
  fi
fi


# Show what we have
echo "resultsDirectory $resultsDirectory"
echo "startIndex: $startIndex"
echo "endIndex: $endIndex"


# Collect list of sides to do (left and/or right)
sides=""
if ( $doLeft ); then
  sides="$sides left"
fi
if ( $doRight ); then
  sides="$sides right"
fi



# Loop over all subjects
for i in `seq $startIndex $endIndex`; do
  # Get subject name
  let subjectIndex=$i-1
  subjectName=${subjectNames[ $subjectIndex ]}


  # Loop over each side
  for side in $sides
  do
    echo "Doing subject $subjectName (index $i) $side"


    # Go into subject left/right directory
    doIt "cd $subjectName/$side"

    # Have a look at the affine registration result
    if ( $doAffineRegistration ); then

      if ( false ); then
        doIt "kvlResample aseg_thresholded.mgz imageDump_coregistered.mgz"
        doIt "kvlViewImage aseg_thresholded_resampled.mgz imageDump_coregistered.mgz"
      else
        doIt "kvlResample nu.mgz imageDump_coregistered.mgz"
        doIt "kvlViewImage nu_resampled.mgz imageDump_coregistered.mgz"
      fi

    fi

    # Have a look at the deformable registration result
    if ( $doDeformableRegistration ); then


      doIt "kvlResample nu.mgz imageDump_coregistered.mgz"
      if ( true ); then
        doIt "kvlResample aseg_thresholded.mgz imageDump_coregistered.mgz"
        doIt "kvlViewImage nu_resampled.mgz aseg_thresholded_resampled.mgz"
      fi

      doIt "kvlViewImage nu_resampled.mgz deformableRegistrationLog/rasterized.mgz"

    fi


    # Have a look at the non-partial volume segmentation result
    if ( $doNonPartialVolumedSegmentation ); then
      doIt "kvlViewImage \
            segmentationWithoutPartialVolumingLog/imageBeingSegmented.mhd \
            segmentationWithoutPartialVolumingLog/colorCodedProbabilies.vtk"
    fi

    # Have a look at the partial volume segmentation result
    if ( $doPartialVolumedSegmentation ); then
      doIt "kvlViewImage \
            segmentationWithPartialVolumingLog/superResolutionImage.mhd \
            segmentationWithPartialVolumingLog/colorCodedProbabilies.vtk"
    fi


    # Move out of directory
    doIt "cd ../.."

  done

done

