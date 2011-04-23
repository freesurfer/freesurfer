#!/bin/bash
#


# Which parts are you going to do
doAffineRegistration=true;
doDeformableRegistration=true;
doNonPartialVolumedSegmentation=true;
doPartialVolumedSegmentation=false;

# K to use
K=0.01;


# Parse input
if [ $# != 4 ]; then
  echo
  echo "Usage: $0 subjectName side inputDirectory outputDirectory"
  echo
  exit -1
fi
subjectName=$1
side=$2
inputDirectory=$3
outputDirectory=$4


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


# Set hard coded location of Atlas3D binaries and atlas
doIt "source kvlSetHardCodedLocations.sh"
echo "atlasDirectory: "$atlasDirectory""

meshFileName="$atlasDirectory/CurrentMeshCollection30.gz"
boundingBoxFileName="$atlasDirectory/imageDump.mgz"


# Depending on the side, retreive the standard FreeSurfer codes for the subfields (cf. FreeSurferColorLUT.txt)
if [ $side = right ]; then
  echo "Doing right side"
  CSF=24
  CerebralWhiteMatter=41
  CerebralCortex=42
  Hippocampus=53
  CA23=500
  CA1=502
  fimbria=503
  presubiculum=504
  hippocampal_fissure=505
  CA4DG=506
  subiculum=507

  compressionLookupTableFileName="$atlasDirectory/compressionLookupTable.txt"

else 
  if [ $side = left ]; then
    echo "Doing left side"
    CSF=24
    CerebralWhiteMatter=2
    CerebralCortex=3
    Hippocampus=17
    CA23=550
    CA1=552
    fimbria=553
    presubiculum=554
    hippocampal_fissure=555
    CA4DG=556
    subiculum=557

    compressionLookupTableFileName="$atlasDirectory/compressionLookupTable_left.txt"

  else 
    echo "$side should be either \"left\" or \"right\""
  fi
fi 



# Show what we have
echo "meshFileName: $meshFileName"
echo "compressionLookupTableFileName: $compressionLookupTableFileName"
echo "boundingBoxFileName: $boundingBoxFileName"
echo "subjectName: $subjectName"
echo "side: $side"
echo "inputDirectory: $inputDirectory"
echo "outputDirectory: $outputDirectory"


# Create a directory to work in, and move into it
workingDirectory="$outputDirectory/$subjectName"
if [ -d $workingDirectory ]; then
  echo "directory $workingDirectory already exists"
else
  doIt "mkdir -p $workingDirectory"
fi
workingDirectory="$workingDirectory/$side"
if [ -d $workingDirectory ]; then
  echo "directory $workingDirectory already exists"
else
  doIt "mkdir -p $workingDirectory"
fi
doIt "cd $workingDirectory"



# Copy aseg.mgz, nu.mgz and talairach.xfm here
doIt "cp -f $inputDirectory/$subjectName/mri/nu.mgz ."
doIt "cp -f $inputDirectory/$subjectName/mri/aseg.mgz ."
doIt "cp -f $inputDirectory/$subjectName/mri/transforms/talairach.xfm ."

# Also copy atlas components here
doIt "cp $boundingBoxFileName ."



##########################################################
# Part I: affine registration
##########################################################

if ( $doAffineRegistration ); then

  if ( false ); then

    # Extract the hippocampus from the aseg segmentation
    doIt "kvlBinaryThresholdImage aseg.mgz $Hippocampus"

    # Coregister the hippocampal prior image affinely with the binary aseg hippocampus
    doIt "kvlRegister imageDump.mgz aseg_thresholded.mgz 12 2"

  else

    if [ $side = left ]; then
      # Apply our right to left atlas" transform (pre-computed)
      doIt "kvlApplyTransform imageDump.mgz \
                              -0.9991   -0.0328    0.0258  144.0010 \
                              -0.0330    0.9994   -0.0086    2.6560 \
                               0.0255    0.0094    0.9996   -3.1939";
      doIt "rm -f imageDump.mgz"
      doIt "mv imageDump_transformed.mgz imageDump.mgz"
    fi

    # Apply "our template"-to-"FS atlas" transform (pre-computed)
    doIt "kvlApplyTransform imageDump.mgz \
                            1.0932    0.0097   -0.0146  -76.4657 \
                           -0.0027    1.0489    0.1291 -121.4683 \
                            0.0153 -0.1434    1.1143  -62.7038";

    # Read the talairach.xfm transform from file
    transformLine=`tail -3 talairach.xfm`;
    tranformLine=${transformLine%;*}; # Separate left part according to ";" symbol

    # Apply the transform
    command=`echo kvlApplyTransform imageDump_transformed.mgz "$tranformLine" 1`;
    doIt "$command"

    # Rename imageDump_transformed_transformed.mgz to imageDump.mgz for historical reasons
    doIt "rm -f imageDump.mgz"
    doIt "mv imageDump_transformed_transformed.mgz imageDump.mgz"

    # Extract the hippocampus from the aseg segmentation
    doIt "kvlBinaryThresholdImage aseg.mgz $Hippocampus"

    # Coregister the hippocampal prior image affinely with the binary aseg hippocampus
    doIt "kvlRegister imageDump.mgz aseg_thresholded.mgz 12 2"

  fi

else
  echo "Skipping affine registration part"
fi



##########################################################
# Part II: deformable registration
##########################################################
deformableRegistrationLogDirectory="deformableRegistrationLog"

if ( $doDeformableRegistration ); then

  # Generate a pseudo intensity image for the aseg
  doIt "kvlBinaryThresholdImage aseg.mgz $Hippocampus $Hippocampus 255 1"

  # Generate a configuration file
  configurationFileName="configurationFileDeformableRegistration.txt"
  doIt "echo \# Configuration for deforming toward aseg                        > $configurationFileName" # Note the > rather than >>!
  doIt "echo logDirectory: $deformableRegistrationLogDirectory                >> $configurationFileName"
  doIt "echo imageFileName: aseg_thresholded.mgz                              >> $configurationFileName"
  doIt "echo meshCollectionFileName: $meshFileName                            >> $configurationFileName"
  doIt "echo K: $K                                                            >> $configurationFileName"
  doIt "echo compressionLookupTableFileName: $compressionLookupTableFileName  >> $configurationFileName"
  doIt "echo boundingFileName: imageDump_coregistered.mgz                     >> $configurationFileName"
  # Force all hippocampal subfields to constitute one class
  doIt "echo sameGaussianParameters: $Hippocampus $CA23 $CA1 $fimbria $presubiculum $hippocampal_fissure $CA4DG $subiculum          >> $configurationFileName"
  # Force all the rest to consitute another class
  doIt "echo sameGaussianParameters: $CSF $CerebralWhiteMatter $CerebralCortex  >> $configurationFileName"
  doIt "echo imageSmoothingSigmas: 0                                          >> $configurationFileName"
  doIt "echo meshSmoothingSigmas: 3                                           >> $configurationFileName"

  # Let the beast go
  doIt "kvlSegmentWithoutGUI $configurationFileName"

else
  echo "Skipping deformable registration part"
fi



##########################################################
# Part III: non-partial volumed segmentation
##########################################################

nonPartialVolumeSegmentationLogDirectory="segmentationWithoutPartialVolumingLog"

# Start with original image as input
inputImageFileName="nu.mgz"

# Possibly do bias field correction as pre-processing step
if ( false ); then

  # Create an image with all ones that describes the bounding box
  doIt "kvlBinaryThresholdImage imageDump_coregistered.mgz 0 255"

  # Resample to get binary mask of all voxels, in native space, inside the cuboid ROI
  doIt "kvlResample imageDump_coregistered_thresholded.mgz nu.mgz 0 0 0"

  # Mask the original image 
  doIt "kvlMaskImage "$inputImageFileName" imageDump_coregistered_thresholded_resampled.mgz"

  # Bias field correct the data
  inputImageFileNameBase=${inputImageFileName%\.*};
  doIt "kvlEMSegment ""$inputImageFileNameBase"_masked.mgz" 3 1"

  # Use the result as input
  inputImageFileName="biasCorrected.mgz"

fi


# Possibly mask out anything far away from aseg hippo segmentation
if ( true ); then

  # Get binary hippo mask from aseg segmentation
  doIt "kvlBinaryThresholdImage aseg.mgz $Hippocampus"

  # Dilate a little bit
  doIt "kvlMathematicalMorphology aseg_thresholded.mgz dilate 3"

  # Mask the original image with the dilated mask
  doIt "kvlMaskImage "$inputImageFileName" aseg_thresholded_dilate_radius3.mgz"

  # Use the result as input
  inputImageFileNameBase=${inputImageFileName%\.*};
  inputImageFileName=""$inputImageFileNameBase"_masked.mgz"
fi
echo "Using $inputImageFileName"





if ( $doNonPartialVolumedSegmentation ); then

  # Upsample the original image
  doIt "kvlUpsample $inputImageFileName 2 2 2"

  # Reconstruct the name of the upsampled image
  inputImageFileNameBase=${inputImageFileName%\.*};
  upsampledImageFileName=""$inputImageFileNameBase"_upsampled_2_2_2.mgz"

  # Generate a configuration file
  configurationFileName="configurationFileSegmentationWithoutPartialVolumingButAfterUpsampling.txt"
  startingMeshCollectionFileName="$deformableRegistrationLogDirectory/warpedOriginalMesh.txt.gz"
  doIt "echo \# Configuration without partial voluming                         > $configurationFileName" # Note the > rather than >>!
  doIt "echo logDirectory: $nonPartialVolumeSegmentationLogDirectory          >> $configurationFileName"
  doIt "echo imageFileName: $upsampledImageFileName                           >> $configurationFileName"
  doIt "echo meshCollectionFileName: $startingMeshCollectionFileName          >> $configurationFileName"
  doIt "echo startingMeshNumber: 0                                            >> $configurationFileName"
  doIt "echo K: $K                                                            >> $configurationFileName"
  doIt "echo compressionLookupTableFileName: $compressionLookupTableFileName  >> $configurationFileName"
  doIt "echo boundingFileName: imageDump_coregistered.mgz                     >> $configurationFileName"
  # doIt "echo logTransformData: 1                                              >> $configurationFileName"
  # doIt "echo biasFieldOrder: 4                                                >> $configurationFileName"
  # Force all WM to have same intensity distributions
  doIt "echo sameGaussianParameters: $CerebralWhiteMatter $fimbria          >> $configurationFileName"
  # Force all GM to have same intensity distributions
  doIt "echo sameGaussianParameters: $CerebralCortex $Hippocampus $CA23 $CA1 $presubiculum $CA4DG $subiculum >> $configurationFileName"
  # Force all CSF to have same intensity distributions
  doIt "echo sameGaussianParameters: $CSF $hippocampal_fissure                >> $configurationFileName"
  doIt "echo imageSmoothingSigmas: 0 0 0                                      >> $configurationFileName"
  doIt "echo meshSmoothingSigmas: 1 0.5 0                                     >> $configurationFileName"

  # Let the beast go
  doIt "kvlSegmentWithoutGUI $configurationFileName"

else
  echo "Skipping non partial volumed segmentation part"
fi


##########################################################
# Part IV: partial volume segmentation
##########################################################

if ( $doPartialVolumedSegmentation ); then

  # Generate a configuration file
  configurationFileName="configurationFileSegmentationWithPartialVoluming.txt"
  partialVolumeSegmentationLogDirectory="segmentationWithPartialVolumingLog"
  startingMeshCollectionFileName="$nonPartialVolumeSegmentationLogDirectory/warpedOriginalMesh.txt.gz"
  doIt "echo \# Configuration with partial voluming                            > $configurationFileName" # Note the > rather than >>!
  doIt "echo logDirectory: $partialVolumeSegmentationLogDirectory             >> $configurationFileName"
  doIt "echo imageFileName: $inputImageFileName                               >> $configurationFileName"
  doIt "echo meshCollectionFileName: $startingMeshCollectionFileName          >> $configurationFileName"
  doIt "echo startingMeshNumber: 0                                            >> $configurationFileName"
  doIt "echo K: $K                                                            >> $configurationFileName"
  doIt "echo compressionLookupTableFileName: $compressionLookupTableFileName  >> $configurationFileName"
  doIt "echo boundingFileName: imageDump_coregistered.mgz                     >> $configurationFileName"
  doIt "echo coregisterToPosteriorProbabilities: 1                            >> $configurationFileName"
  doIt "echo partialVolumeUpsamplingFactors: 2 2 2                            >> $configurationFileName"
  # Force all WM to have same intensity distributions
  doIt "echo sameGaussianParameters: $CerebralWhiteMatter $fimbria          >> $configurationFileName"
  # Force all GM to have same intensity distributions
  doIt "echo sameGaussianParameters: $CerebralCortex $Hippocampus $CA23 $CA1 $presubiculum $CA4DG $subiculum >> $configurationFileName"
  # Force all CSF to have same intensity distributions
  doIt "echo sameGaussianParameters: $CSF $hippocampal_fissure                >> $configurationFileName"
  # doIt "echo imageSmoothingSigmas: 0 0 0                                      >> $configurationFileName"
  # doIt "echo meshSmoothingSigmas: 1 0.5 0                                     >> $configurationFileName"

  # Let the beast go
  doIt "kvlSegmentWithoutGUI $configurationFileName"

else
  echo "Skipping the partial volume segmentation part"
fi



