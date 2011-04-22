#!/bin/bash

#
# ./doAllSubjects.sh /autofs/space/madrc_007/users/jab/FS4Koen/ alex 1 4
#


# Which parts are you going to do
doLeft=true;
doRight=true;



# Parse input
if [ $# -lt 2 ]; then
  echo
  echo "Usage: $0 inputDirectory outputDirectory [ startIndex=1 endIndex=end ]"
  echo
  echo "    e.g.:"
  echo
  echo "       $0 /homes/13/koen/erdosSpace/data/alexBecker/FS4Koen/ /homes/13/koen/erdosSpaceDisk4/data/alex 1 4"
  echo
  exit -1
fi
inputDirectory=$1
outputDirectory=$2

startIndex=1
if [ $# -ge 3 ]; then
  startIndex=$3
fi

subjectNames=(`ls $inputDirectory`)
if [ $? != 0 ]; then
  echo "failed to list file names in $inputDirectory"
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
echo "inputDirectory: $inputDirectory"
echo "outputDirectory $outputDirectory"
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


# Create the output directory if it doesn't yet exist
if [ -d $outputDirectory ]; then
  echo "directory $outputDirectory already exists"
else
  doIt "mkdir $outputDirectory"
fi


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
    echo "Doing subject $subjectName $side"

    # Construct output file name (catches both cout and cerr)
    logFileName="$outputDirectory/logOfSubject_"$subjectName"_"$side".txt"
    echo "logFileName: $logFileName"

    # Let the beast go
    doIt "doSubject.sh $subjectName $side $inputDirectory $outputDirectory &> $logFileName"
  done

done

