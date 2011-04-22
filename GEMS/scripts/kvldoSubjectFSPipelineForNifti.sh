#!/bin/bash
#
# Let's do one subject
#


# Parse input
if [ $# != 3 ]; then
  echo
  echo "Usage: $0 subjectName originalDataDirectory freeSurferResultsDirectory"
  echo
  echo
  echo "    e.g.:"
  echo
  echo "       $0 KVL_A05_anat \
                  /homes/13/koen/erdosSpaceDisk4/data/matthewMalterCohen \
                  /homes/13/koen/erdosSpaceDisk4/data/matthewMalterCohenFreeSurferResults"
  echo
  exit -1
fi
subjectName=$1
originalDataDirectory=$2
freeSurferResultsDirectory=$3


# Show what we have
echo "subjectName: $subjectName"
echo "originalDataDirectory: $originalDataDirectory"
echo "freeSurferResultsDirectory: $freeSurferResultsDirectory"


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
workingDirectory="$freeSurferResultsDirectory/"$subjectName"/mri/orig"
if [ -d $workingDirectory ]; then
        echo "directory $workingDirectory already exists"
else
        doIt "mkdir -p $workingDirectory"
fi
doIt "cd $workingDirectory"


# Convert Nifti to MGZ
niftiFileName="$originalDataDirectory/"$subjectName".nii"
doIt "mri_convert $niftiFileName 001.mgz"

# Start FreeSurfer pipeline
export SUBJECTS_DIR=$freeSurferResultsDirectory
doIt "recon-all -all -s $subjectName"


