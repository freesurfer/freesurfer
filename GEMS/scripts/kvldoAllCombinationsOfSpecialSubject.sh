#!/bin/bash

#
# ./doAllSubjects.sh /autofs/space/madrc_007/users/jab/FS4Koen/ alex 1 4
#

# Parse input
if [ $# -lt 3 ]; then
  echo
  echo "Usage: $0 inputDirectory outputDirectory subjectName"
  echo
  echo "    e.g.:"
  echo
  echo "       $0 /homes/13/koen/erdosSpace/data/alexBecker/FS4Koen/ /homes/13/koen/erdosSpaceDisk4/data/special 341-1001"
  echo
  exit -1
fi
inputDirectory=$1
outputDirectory=$2
subjectName=$3


# Show what we have
echo "inputDirectory: $inputDirectory"
echo "outputDirectory $outputDirectory"
echo "subjectName: $subjectName"


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


# Loop over all combinations of K and meshSmoothingSigma
startK=0.01
endK=1.0
numberOfKSteps=8

startSigma=0
endSigma=3
numberOfSigmaSteps=6


logStartK=$(echo "l($startK)" | bc -l)
logEndK=$(echo "l($endK)" | bc -l)
intervalLogK=$(echo "( $logEndK - $logStartK ) / ( $numberOfKSteps - 1 )" | bc -l)

intervalSigma=$(echo "( $endSigma - $startSigma ) / ( $numberOfSigmaSteps - 1 )" | bc -l)

for logK in `seq -f %f $logStartK $intervalLogK $logEndK`; do
  # Compute K
  K=$(echo "scale=2;e($logK)" | bc -l)
  #echo $K

  for sigma in `seq $startSigma $intervalSigma $endSigma`; do
    sigma=$(echo "scale=2;e(l($sigma))" | bc -l)
    #echo $sigma

    # Display what we're doing
    echo "Doing combination K="$K" sigma="$sigma""

    # Construct output file name (catches both cout and cerr)
    logFileName="$outputDirectory/logOfSubject_"$subjectName"__K_"$K"__meshSmoothingSigma_"$sigma".txt"
    echo "logFileName: $logFileName"

    # Let the beast go
    doIt "./doSpecialSubject.sh $subjectName $inputDirectory $outputDirectory $K $sigma &> $logFileName"

  done

done
