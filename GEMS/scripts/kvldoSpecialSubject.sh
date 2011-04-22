#!/bin/bash
#
# Let's do one subject
#

#
# ./doSubject.sh 129-1006 /autofs/space/madrc_007/users/jab/FS4Koen/ alex 0.01 0.5
#


# Which parts are you going to do
doAffineRegistration=true;
doNonPartialVolumedSegmentation=true;
doPartialVolumedSegmentation=false;


# Hard coded location of Atlas3D binaries
PATH=$PATH:/homes/13/koen/erdosSpace/software/Atlas3D-Release-CompiledOnOct/bin

# Hard coded locations of atlas
atlasDirectory=/homes/13/koen/erdosSpace/software/Atlas3D-Release-CompiledOnOct/bin/atlas
meshFileName="$atlasDirectory/CurrentMeshCollection30.gz"
compressionLookupTableFileName="$atlasDirectory/compressionLookupTable.txt"
boundingBoxFileName="$atlasDirectory/imageDump.mgz"


# Parse input
if [ $# != 5 ]; then
  echo
  echo "Usage: $0 subjectName inputDirectory outputDirectory K meshSmoothingSigma"
  echo
  echo
  echo "    e.g.:"
  echo
  echo "       $0 341-1001 /homes/13/koen/erdosSpace/data/alexBecker/FS4Koen/ /homes/13/koen/erdosSpaceDisk4/data/special 0.02 0.5"
  echo
  exit -1
fi
subjectName=$1
inputDirectory=$2
outputDirectory=$3
K=$4
meshSmoothingSigma=$5


# Show what we have
echo "meshFileName: $meshFileName"
echo "compressionLookupTableFileName: $compressionLookupTableFileName"
echo "boundingBoxFileName: $boundingBoxFileName"
echo "subjectName: $subjectName"
echo "inputDirectory: $inputDirectory"
echo "outputDirectory: $outputDirectory"
echo "K $K"
echo "meshSmoothingSigma: $meshSmoothingSigma"


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


# Create a directory to work in, and move into it
workingDirectory="$outputDirectory/"$subjectName"__K_"$K"__meshSmoothingSigma_"$meshSmoothingSigma""
if [ -d $workingDirectory ]; then
        echo "directory $workingDirectory already exists"
else
        doIt "mkdir $workingDirectory"
fi
doIt "cd $workingDirectory"


# Copy aseg.mgz and nu.mgz here
doIt "cp $inputDirectory/$subjectName/mri/nu.mgz ."
doIt "cp $inputDirectory/$subjectName/mri/aseg.mgz ."

# Also copy atlas components here
doIt "cp $boundingBoxFileName ."



##########################################################
# Part I: affine registration
##########################################################

if ( $doAffineRegistration ); then

  # Extract the right hippocampus from the aseg segmentation (label 53)
  doIt "kvlBinaryThresholdImage aseg.mgz 53"

  # Coregister the hippocampal prior image affinely with the binary aseg hippocampus
  doIt "kvlRegister imageDump.mgz aseg_thresholded.mgz 12 2"

else
  echo "Skipping affine registration part"
fi


##########################################################
# Part II: non-partial volumed segmentation
##########################################################

nonPartialVolumeSegmentationLogDirectory="segmentationWithoutPartialVolumingLog"

if ( $doNonPartialVolumedSegmentation ); then

  # Upsample the original image
  doIt "kvlUpsample nu.mgz 2 2 2"

  # Generate a configuration file
  configurationFileName="configurationFileSegmentationWithoutPartialVolumingButAfterUpsampling.txt"
  doIt "echo \# Configuration without partial voluming                         > $configurationFileName" # Note the > rather than >>!
  doIt "echo logDirectory: $nonPartialVolumeSegmentationLogDirectory          >> $configurationFileName"
  doIt "echo imageFileName: nu_upsampled_2_2_2.mgz                            >> $configurationFileName"
  doIt "echo meshCollectionFileName: $meshFileName                            >> $configurationFileName"
  # doIt "echo K: 0.01                                                          >> $configurationFileName"
  doIt "echo K: $K                                                          >> $configurationFileName"
  doIt "echo compressionLookupTableFileName: $compressionLookupTableFileName  >> $configurationFileName"
  doIt "echo boundingFileName: imageDump_coregistered.mgz                     >> $configurationFileName"
  # Force all WM to have same intensity distributions
  doIt "echo sameGaussianParameters: 41 503                                   >> $configurationFileName"
  # Force all GM to have same intensity distributions
  doIt "echo sameGaussianParameters: 42 53 500 502 504 506 507                >> $configurationFileName"
  # Force all CSF to have same intensity distributions
  doIt "echo sameGaussianParameters: 24 505                                   >> $configurationFileName"
  # doIt "echo imageSmoothingSigmas: 0 0 0                                      >> $configurationFileName"
  # doIt "echo meshSmoothingSigmas: 1 0.5 0                                     >> $configurationFileName"
  doIt "echo imageSmoothingSigmas: 0                                          >> $configurationFileName"
  doIt "echo meshSmoothingSigmas: $meshSmoothingSigma:                         >> $configurationFileName"


  # Let the beast go
  doIt "kvlSegmentWithoutGUI $configurationFileName"

else
  echo "Skipping non partial volumed segmentation part"
fi


##########################################################
# Part III: partial volume segmentation
##########################################################

if ( $doPartialVolumedSegmentation ); then

  # Generate a configuration file
  configurationFileName="configurationFileSegmentationWithPartialVoluming.txt"
  partialVolumeSegmentationLogDirectory="segmentationWithPartialVolumingLog"
  startingMeshCollectionFileName="$nonPartialVolumeSegmentationLogDirectory/warpedOriginalMesh.txt.gz"
  doIt "echo \# Configuration with partial voluming                            > $configurationFileName" # Note the > rather than >>!
  doIt "echo logDirectory: $partialVolumeSegmentationLogDirectory             >> $configurationFileName"
  doIt "echo imageFileName: nu.mgz                                            >> $configurationFileName"
  doIt "echo meshCollectionFileName: $startingMeshCollectionFileName          >> $configurationFileName"
  doIt "echo startingMeshNumber: 0                                            >> $configurationFileName"
  doIt "echo K: 0.01                                                          >> $configurationFileName"
  doIt "echo compressionLookupTableFileName: $compressionLookupTableFileName  >> $configurationFileName"
  doIt "echo boundingFileName: imageDump_coregistered.mgz                     >> $configurationFileName"
  doIt "echo coregisterToPosteriorProbabilities: 1                            >> $configurationFileName"
  doIt "echo partialVolumeUpsamplingFactors: 2 2 2                            >> $configurationFileName"
  # Force all WM to have same intensity distributions
  doIt "echo sameGaussianParameters: 41 503                                   >> $configurationFileName"
  # Force all GM to have same intensity distributions
  doIt "echo sameGaussianParameters: 42 53 500 502 504 506 507                >> $configurationFileName"
  # Force all CSF to have same intensity distributions
  doIt "echo sameGaussianParameters: 24 505                                   >> $configurationFileName"
  # doIt "echo imageSmoothingSigmas: 0 0 0                                      >> $configurationFileName"
  # doIt "echo meshSmoothingSigmas: 1 0.5 0                                     >> $configurationFileName"

  # Let the beast go
  doIt "kvlSegmentWithoutGUI $configurationFileName"

else
  echo "Skipping the partial volume segmentation part"
fi



