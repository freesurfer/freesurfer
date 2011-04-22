#!/bin/bash

# Which parts are you going to collect
collectLeft=true;
collectRight=true;
collectDeformableRegistration=true;
collectNonPartialVolumedSegmentation=true;
collectPartialVolumedSegmentation=false;


# Parse input
if [ $# -lt 1 ]; then
  echo
  echo "Usage: $0 outputDirectory"
  echo
  echo "    e.g.:"
  echo
  echo "       $0 ~/erdosSpaceDisk4/data/alex"
  echo
  exit -1
fi
outputDirectory=$1


# Show what we have
echo "outputDirectory $outputDirectory"


# Define a function that will show what we're doing, and
# exits if something goes wrong.
function doIt {

  command="$1"

  echo "$command"
  eval "$command"

  if [ $? != 0 ]; then
    echo "failed to do "$command""
  exit -1
fi

}



# Set hard coded location of Atlas3D binaries and atlas
doIt "source kvlSetHardCodedLocations.sh"
echo "atlasDirectory: "$atlasDirectory""


# Collect list of sides to do (left and/or right)
sides=""
if ( $collectLeft ); then
  sides="$sides left"
fi
if ( $collectRight ); then
  sides="$sides right"
fi



# Go to the output directory
doIt "cd $outputDirectory"

# List all subdirectories
subjectNames=(`ls -d */`) # list only directories
if [ $? != 0 ]; then
  echo "failed to list file names in $outputDirectory"
  exit -1
fi
numberOfSubjects=${#subjectNames[*]}


# Go inside each non-partial and partial volume segmentation subdirectory
# and compose a color coded segmentation
for (( i = 0; i < $numberOfSubjects; i++ )); do
  # Get subject name
  subjectName=${subjectNames[ $i ]}

  # Loop over each side
  for side in $sides; do
    echo "Collecting subject $subjectName $side"

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

      rasterizeLeftRightFlag=1;

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

      rasterizeLeftRightFlag=2;

    fi 


    # Get inside the subject's right/left directory
    doIt "cd $subjectName/$side"

    if ( $collectDeformableRegistration ); then

      # Get inside the deformable registration
      doIt "cd deformableRegistrationLog"

      # Rasterize the prior for global hippocampus
      doIt "kvlMapLabelsOfAtlasMesh warpedOriginalMesh.txt.gz $compressionLookupTableFileName $rasterizeLeftRightFlag"
      doIt "kvlRasterizeAtlasMesh warpedOriginalMesh.txt.gz_mapped.txt.gz \
            ../imageDump_coregistered.mgz 1 0"

      # Go back to the subject's main directory
      doIt "cd .."

    else
      echo "Skipping deformable registration part"
    fi


    if ( $collectNonPartialVolumedSegmentation ); then

      # Get inside the non-partial volume subdir
      doIt "cd segmentationWithoutPartialVolumingLog"

      # Color code the segmentation
      doIt "kvlColorCodeProbabilityImages $compressionLookupTableFileName \
            posterior_"$Hippocampus".mgz \
            posterior_"$presubiculum".mgz \
            posterior_"$CA1".mgz \
            posterior_"$CA23".mgz \
            posterior_"$fimbria".mgz \
            posterior_"$subiculum".mgz \
            posterior_"$CA4DG".mgz \
            posterior_"$hippocampal_fissure".mgz"

      # Go back to the subject's main directory
      doIt "cd .."

    else
      echo "Skipping non-partial volume segmentation part"
    fi

    if ( $collectPartialVolumedSegmentation ); then

      # Get inside the partial volume subdir
      doIt "cd segmentationWithPartialVolumingLog"

      # Color code the segmentation
      doIt "kvlColorCodeProbabilityImages $compressionLookupTableFileName \
            superResolutionPosterior_"$Hippocampus".mgz \
            superResolutionPosterior_"$presubiculum".mgz \
            superResolutionPosterior_"$CA1".mgz \
            superResolutionPosterior_"$CA23".mgz \
            superResolutionPosterior_"$fimbria".mgz \
            superResolutionPosterior_"$subiculum".mgz \
            superResolutionPosterior_"$CA4DG".mgz \
            superResolutionPosterior_"$hippocampal_fissure".mgz"

      # Go back to the subject's main directory
      doIt "cd .."

    else
      echo "Skipping partial volume segmentation part"
    fi


    # Go back to main directory
    doIt "cd ../.."

  done

done



# Collect results
if $collectLeft && $collectRight; then
  wildCardString="{left,right}"
else

  if ( $collectLeft ); then
    wildCardString="left"
  else
    wildCardString="right"
  fi

fi
# echo $wildCardString


command="tar cvfz results.tgz \
      */"$wildCardString"/aseg_thresholded.mgz \
      */"$wildCardString"/nu.mgz \
      */"$wildCardString"/imageDump_coregistered.mgz"
if ( $collectDeformableRegistration ); then
  command=""$command" \
           */"$wildCardString"/deformableRegistrationLog/rasterized.mgz"
fi
if ( $collectNonPartialVolumedSegmentation ); then
  command=""$command" \
           */"$wildCardString"/segmentationWithoutPartialVolumingLog/colorCodedProbabilies.vtk \
           */"$wildCardString"/segmentationWithoutPartialVolumingLog/imageBeingSegmented.*"
fi
if ( $collectPartialVolumedSegmentation ); then
  command=""$command" \
           */"$wildCardString"/segmentationWithPartialVolumingLog/colorCodedProbabilies.vtk \
           */"$wildCardString"/segmentationWithPartialVolumingLog/superResolutionImage.*"
fi

doIt "$command"


# Show where we've written results
echo
echo "results collected in $outputDirectory/results.tgz"
echo


