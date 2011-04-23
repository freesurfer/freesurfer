#!/bin/bash


# Parse input
if [ $# -lt 2 ]; then
  echo
  echo "Usage: $0 originalDataDirectory freeSurferResultsDirectory [ startIndex=1 endIndex=end ]"
  echo
  echo "    e.g.:"
  echo
  echo "       $0 /homes/13/koen/erdosSpaceDisk4/data/matthewMalterCohen \
                  /homes/13/koen/erdosSpaceDisk4/data/matthewMalterCohenFreeSurferResults 1 4"
  echo
  exit -1
fi
originalDataDirectory=$1
freeSurferResultsDirectory=$2

startIndex=1
if [ $# -ge 3 ]; then
  startIndex=$3
fi

subjectNames=(`ls $originalDataDirectory | grep .nii`)
if [ $? != 0 ]; then
  echo "failed to list Nifti file names in $originalDataDirectory"
  exit -1
fi
numberOfSubjects=${#subjectNames[*]}
#echo "numberOfSubjects: $numberOfSubjects"

endIndex=$numberOfSubjects
if [ $# -ge 4 ]; then
  endIndex=$4
  if [ $endIndex -gt $numberOfSubjects ]; then
    echo "endIndex ($endIndex) can't be greater than the number of subjects ($numberOfSubjects)"
    exit -1
  fi
fi


# Show what we have
echo "originalDataDirectory: $originalDataDirectory"
echo "freeSurferResultsDirectory $freeSurferResultsDirectory"
echo "startIndex: $startIndex"
echo "endIndex: $endIndex"


# Define a function that will show what we're doing, and
# exits if something goes wrong. It also detects if we're
# running on seychelles, and if so will submit the input
# argument as a job. Otherwise the argument is just run locally
function doIt {

  command=$1
  if [ `hostname` = seychelles ]; then
    command="pbsubmit -l nodes=1:opteron -c \""$command"\""
  fi

  echo $command
  eval $command

  if [ $? != 0 ]; then
    echo "failed to do $command"
  exit -1
fi

}


# Loop over all subjects
for i in `seq $startIndex $endIndex`; do
  # Get subject name
  let subjectIndex=$i-1
  subjectName=${subjectNames[ $subjectIndex ]}
  subjectName=${subjectName%.nii} # Remove trailing .nii
  echo "Doing subject $subjectName"

  # Construct output file name (catches both cout and cerr)
  logFileName="$freeSurferResultsDirectory/logOfSubject_"$subjectName".txt"
  echo "logFileName: $logFileName"

  # Let the beast go
  doIt "./doSubjectFSPipelineForNifti.sh $subjectName $originalDataDirectory $freeSurferResultsDirectory &> $logFileName"

done

