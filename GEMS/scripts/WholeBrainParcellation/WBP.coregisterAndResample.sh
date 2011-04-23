#!/bin/bash
#
# Set up all the coregistration and resampling stuff we need
#
# If this script failed and was stopped, do 
# 
#  rm *.mgz; rm FreeSurferColorLUT.txt; rm subject*/*coregistered*; rm subject*/FreeSurferColorLUT.txt
# 
# to clean up the mess and start with a clean plate
#


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


# Make sure we have FreeSurfer's MNI template
doIt "ln -s /usr/local/freesurfer/dev/average/mni305.cor.mgz ."
doIt "ln -s /usr/local/freesurfer/dev/average/mni305.mask.cor.mgz ."

# Also copy FreeSurfer's lookup table
doIt "cp /homes/13/koen/software/Atlas3D/FreeSurferColorLUT.txt ."

# Mask out the brain in the MNI template
doIt "kvlMaskImage mni305.cor.mgz mni305.mask.cor.mgz" 
doIt "kvlViewImage mni305_masked.mgz"

# Auto crop the brain masked MNI template, adding a 10 voxel wide
# border. This will be the image grid on which our mesh-based 
# atlases are built, while still being in the correct MNI space
# through the image-to-world transform in the masked and bordered
# template
doIt "kvlAutoCrop mni305_masked.mgz 10"
doIt "kvlViewImage mni305_masked_autoCropped.mgz"


# Now visit each subject: list all subdirectories and do work inside each one
subjectNames=(`ls -d */`) # list only directories
if [ $? != 0 ]; then
  echo "failed to list file names in $PWD"
  exit -1
fi
numberOfSubjects=${#subjectNames[*]}
for (( i = 0; i < $numberOfSubjects; i++ )); do
  # Get subject name
  subjectName=${subjectNames[ $i ]}

  # Get inside the subject's directory
  doIt "cd $subjectName"

  # Generate a skull stripping mask using mathematical morphology
  # closing operation on the manual segmenations (which don't 
  # include cortical CSF). We might want to consider eroding with a
  # smaller kernel than the original dilation kernel in order to really
  # include additional CSF (?)
  doIt "kvlBinaryThresholdImage seg_edited.mgz 1 10000"
  doIt "kvlMathematicalMorphology seg_edited_thresholded.mgz dilate 10"
  doIt "kvlMathematicalMorphology seg_edited_thresholded_dilate_radius10.mgz \
        erode 10"
  doIt "ln -s seg_edited_thresholded_dilate_radius10_erode_radius10.mgz mask.mgz"

  # Inspect the mask
  doIt "kvlViewImage orig.mgz mask.mgz"

  # Apply the mask to the data (i.e., skull-strip)
  doIt "kvlMaskImage orig.mgz mask.mgz"

  # Inspect the result
  doIt "kvlThresholdImage orig_masked.mgz"

  # Coregister the skull-stripped image with our preprocessed MNI 
  # template, starting by aligning the momenta and applying the
  # same transform to the manual labels
  doIt "kvlRegister orig_masked.mgz ../mni305_masked_autoCropped.mgz \
        12 2 20 1 seg_edited.mgz"

  # Resample both the skull-stripped image and the manual segs to
  # the image grid of our preprocessed MNI template. Use nearest-
  # neighbor interpolation.
  doIt "kvlResample orig_masked_coregistered.mgz \
        ../mni305_masked_autoCropped.mgz 0"
  doIt "kvlResample seg_edited_coregistered.mgz \
        ../mni305_masked_autoCropped.mgz 0"

  # Visually inspect the registration result
  doIt "kvlViewImage ../mni305_masked_autoCropped.mgz \
        orig_masked_coregistered_resampled.mgz" 

  # Also make sure our labels our transformed correctly as well
  doIt "cp ../FreeSurferColorLUT.txt ."
  doIt "kvlCompressImage seg_edited_coregistered_resampled.mgz"
  doIt "kvlViewImage orig_masked_coregistered_resampled.mgz \
        compressed.mhd compressionLookupTable.txt"

  # Now add artificially some cortical CSF to the manual 
  # segmentations, using the voxels that lie in our mask but
  # have not been manually labeled
  doIt "kvlAddCSFToBucknerData seg_edited_coregistered_resampled.mgz \
   orig_masked_coregistered_resampled.mgz"

  # Look at the final result
  doIt "kvlCompressImage BucknerSegmentationWithAddedCSF.mgz"
  doIt "kvlViewImage orig_masked_coregistered_resampled.mgz \
        compressed.mhd compressionLookupTable.txt"

  # Go back to main directory
  doIt "cd .."

done






