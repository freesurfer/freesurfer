#!/bin/bash
#
# Extracts and then co-registers, using affine registration, a certain 
# structure from ASEG segmentations.
#



# Parse input
if [ $# -lt 2 ]; then
  echo
  echo "Usage: $0 inputImage referenceImage [ showResult=false label=54 ] "
  echo
  exit -1
fi
inputImage=$1
referenceImage=$2
showResult=false
if [ $# -ge 3 ]; then
  showResult=$3
fi
label=54;
if [ $# -ge 4 ]; then
  label=$4
fi

# Show what we have
echo "inputImage: $inputImage"
echo "referenceImage: $referenceImage"
echo "showResult: $showResult"
echo "label: $label"


# Define a function that will show what we're doing, and
# exits if something goes wrong
function doIt {
  echo $1
  eval $1
  if [ $? != 0 ]; then
    echo "failed to do $1"
  exit -1
fi

}


# Copy input asegs
doIt "cp $inputImage inputImage_aseg.mgz";
doIt "cp $referenceImage referenceImage_aseg.mgz";

# Extract the right structure from the aseg segmentations
doIt "kvlBinaryThresholdImage inputImage_aseg.mgz $label";
doIt "kvlBinaryThresholdImage referenceImage_aseg.mgz $label";

# Crop the binary input image
doIt "kvlAutoCrop inputImage_aseg_thresholded.mgz 10";

# Coregister the two binary volumes
doIt "kvlRegister \
      inputImage_aseg_thresholded_autoCropped.mgz \
      referenceImage_aseg_thresholded.mgz \
      12 2 20 1 \
      inputImage_aseg_thresholded.mgz";

# Resample the input volume
doIt "kvlResample inputImage_aseg_thresholded_coregistered.mgz referenceImage_aseg_thresholded.mgz";

# Rename
doIt "mv inputImage_aseg_thresholded_coregistered_resampled.mgz inputResampledBinaryVolume.mgz";
doIt "mv referenceImage_aseg_thresholded.mgz referenceBinaryVolume.mgz";

# Clean up
doIt "rm inputImage_aseg.mgz referenceImage_aseg.mgz inputImage_aseg_thresholded.mgz \
         inputImage_aseg_thresholded_autoCropped.mgz inputImage_aseg_thresholded_coregistered.mgz \
         inputImage_aseg_thresholded_autoCropped_coregistered.mgz";

# Show 
if ( $showResult ); then

  doIt "kvlViewImage referenceBinaryVolume.mgz inputResampledBinaryVolume.mgz";

fi





