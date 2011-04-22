#!/bin/bash
#
# "Usage: $0 [ resultsDirectory startIndex=1 endIndex=end ]"
#



# Which parts are you going to quantify
doLeft=true;
doRight=true;
doNonPartialVolumedSegmentation=true;
doPartialVolumedSegmentation=false;


# Define a function that will show what we're doing, and
# exits if something goes wrong.
function doIt {

  command="$1"

  echo "$command"
  eval "$command"

  if [ $? != 0 ]; then
    echo "failed to do $command"
    exit -1
  fi
}


# Set hard coded location of Atlas3D binaries and atlas
source kvlSetHardCodedLocations.sh


# Print help if needed
if [ $# = 1 ]; then
  if [ $1 = "--help" ]; then
    fsPrintHelp $atlasDirectory/kvlQuantifyHippocampalSubfieldSegmentations.sh.help.xml
    exit 0
  fi
fi


# See what type of directory structure this is: the one created by
# kvlSegmentHippocampalSubfields.sh, which has its results in
# something like 
#   $resultsDirectory/bert/left/segmentationWithoutPartialVolumingLog,
# or the FreeSurfer one where those results have been copied inot
#   $SUBJECTS_DIR/bert/mri

resultsDirectory=$SUBJECTS_DIR
useOriginalStructure=false
if [ $# -ge 1 ]; then
  resultsDirectory=$1
  useOriginalStructure=true
fi


# Show what we have
echo "resultsDirectory $resultsDirectory"


# Go to the output directory
doIt "cd $resultsDirectory"

# List all subdirectories
subjectNames=(`ls -d */`) # list only directories
if [ $? != 0 ]; then
  echo "failed to list file names in $resultsDirectory"
  exit -1
fi
numberOfSubjects=${#subjectNames[*]}


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


# Go inside each subject's subdirectory and calculate volume statistics
# there for non-partial volume and/or partial volume segmentations
#for i in `seq $startIndex $endIndex`; do
for i in `eval echo {$startIndex..$endIndex}`; do
  # Get subject name
  let subjectIndex=$i-1
  subjectName=${subjectNames[ $subjectIndex ]}

  # Go into subject's main directory
  doIt "cd $subjectName"

  # Loop over each side
  for side in $sides; do
    echo "Quantifying subject $subjectName $side"

    # Depending on the side, retreive the standard FreeSurfer names for the subfields (cf. FreeSurferColorLUT.txt)
    if [ $side = right ]; then
      echo "Doing right side"
      Hippocampus="Right-Hippocampus"
      CA23="right_CA2-3"
      CA1="right_CA1"
      fimbria="right_fimbria"
      presubiculum="right_presubiculum"
      hippocampal_fissure="right_hippocampal_fissure"
      CA4DG="right_CA4-DG"
      subiculum="right_subiculum"

      compressionLookupTableFileName="$atlasDirectory/compressionLookupTable.txt"

    else 
      echo "Doing left side"
      Hippocampus="Left-Hippocampus"
      CA23="left_CA2-3"
      CA1="left_CA1"
      fimbria="left_fimbria"
      presubiculum="left_presubiculum"
      hippocampal_fissure="left_hippocampal_fissure"
      CA4DG="left_CA4-DG"
      subiculum="left_subiculum"

      compressionLookupTableFileName="$atlasDirectory/compressionLookupTable_left.txt"

    fi 

    # Get inside the right/left directory
    if ( $useOriginalStructure ); then
      doIt "cd $side"
    else
      doIt "cd mri"
    fi

    # Collect statistics
    if ( $doNonPartialVolumedSegmentation ); then

      if ( $useOriginalStructure ); then
        # Get inside the non-partial volume subdir
        doIt "cd segmentationWithoutPartialVolumingLog"
      fi

      # Quantify
      doIt "kvlQuantifyPosteriorProbabilityImages $compressionLookupTableFileName \
            posterior_"$Hippocampus".mgz \
            posterior_"$presubiculum".mgz \
            posterior_"$CA1".mgz \
            posterior_"$CA23".mgz \
            posterior_"$fimbria".mgz \
            posterior_"$subiculum".mgz \
            posterior_"$CA4DG".mgz \
            posterior_"$hippocampal_fissure".mgz \
            > volumeStats_$side.txt"

      if ( $useOriginalStructure ); then
        # Get out of this subdir
        doIt "cd .."
      fi

    fi

    if ( $doPartialVolumedSegmentation ); then

      if ( $useOriginalStructure ); then
        # Get inside the non-partial volume subdir
        doIt "cd segmentationWithPartialVolumingLog"
      fi

      # Quantify
      doIt "kvlQuantifyPosteriorProbabilityImages $compressionLookupTableFileName \
            posterior_"$Hippocampus".mgz \
            posterior_"$presubiculum".mgz \
            posterior_"$CA1".mgz \
            posterior_"$CA23".mgz \
            posterior_"$fimbria".mgz \
            posterior_"$subiculum".mgz \
            posterior_"$CA4DG".mgz \
            posterior_"$hippocampal_fissure".mgz \
            > volumeStats_$side.txt"

      if ( $useOriginalStructure ); then
        # Get out of this subdir
        doIt "cd .."
      fi
    fi

    # Go back into subject's main directory
    doIt "cd .."

  done # End loop over left/right side

  # Go back to main directory
  doIt "cd .."

done # End loop over all subjects



# Now collect the results from all the files

# Loop over each side
for side in $sides; do
  echo "Collecting all results for side: $side"

  # Construct correct results file name
  if [ $side = right ]; then

    if ( $doNonPartialVolumedSegmentation ); then
      allStatsFileName="nonPartialVolumeStatsRight.txt"
    else
      allStatsFileName="partialVolumeStatsRight.txt"
    fi

  else 

    if ( $doNonPartialVolumedSegmentation ); then
      allStatsFileName="nonPartialVolumeStatsLeft.txt"
    else
      allStatsFileName="partialVolumeStatsLeft.txt"
    fi

  fi 



  #for i in `seq $startIndex $endIndex`; do
  for i in `eval echo {$startIndex..$endIndex}`; do
    # Get subject name
    let subjectIndex=$i-1
    subjectName=${subjectNames[ $subjectIndex ]}


    # Parse the result
    if ( $useOriginalStructure ); then
 
      if ( $doNonPartialVolumedSegmentation ); then
        statFileName="$subjectName/$side/segmentationWithoutPartialVolumingLog/volumeStats_$side.txt"
      else
        statFileName="$subjectName/$side/segmentationWithPartialVolumingLog/volumeStats_$side.txt"
      fi

    else
      statFileName="$subjectName/mri/volumeStats_$side.txt"
    fi

    leftPartsString="volume_in_number_of_voxels"
    rightPartsString=${subjectName//\//}  # Substring replacement syntax ${string//substring/replacement}
    while read line; do
      #echo this line is: $line;

      # Separate line according to ":" symbol
      leftPart=${line%:*}
      rightPart=${line#*:}
      #echo "         leftPart: $leftPart"
      #echo "         rightPart: $rightPart"

      #
      if [ $leftPart != "volumeInVoxels" ]; then
        leftPartsString="$leftPartsString   $leftPart"
        rightPartsString="$rightPartsString   $rightPart"
      fi

    done < $statFileName


    # If this is the first subject, start the file and write out the 
    # label names in the first line
    if [ $i -eq $startIndex ]; then
      doIt "echo \"$leftPartsString\"   > $allStatsFileName"
    fi

    # Append a line with the numerical results for this subject to the file
    doIt "echo \"$rightPartsString\"    >> $allStatsFileName"

  done


  # # Show where we've written results
  echo
  echo "results collected in $allStatsFileName"
  echo


done # End loop over left/right side





