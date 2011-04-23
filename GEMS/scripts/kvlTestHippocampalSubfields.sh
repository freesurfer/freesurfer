#!/bin/bash
#


# Parse input
if [ $# != 1 ]; then
  echo
  echo "Usage: $0 inputDirectory"
  echo
  echo "This script will segment the hippocampal subfields for test subject 'bert',"
  echo "calculate the volumes of the subfields, and compare those to previously"
  echo "computed reference values. The 'inputDirectory' is the parent directory of"
  echo "the 'bert' directory with FreeSurfer's standard pipeline segmentation results,"
  echo "i.e., typically \$FREESURFER_HOME/subjects"
  echo
  exit -1
fi
inputDirectory=$1



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

# Make sure we know where to find the reference tables with results
# we're going to compare ourselves to
doIt "source kvlSetHardCodedLocations.sh"
echo "atlasDirectory: "$atlasDirectory""


# Create a directory holding the test results and go inside it
resultDirectory="/tmp/hippocampalSubfieldTestDirectory"
doIt "mkdir -p $resultDirectory"
doIt "cd $resultDirectory"

# Segment both the left and right side of the "bert" test subject
doIt "kvlSegmentHippocampalSubfields.sh bert left $inputDirectory $resultDirectory"
doIt "kvlSegmentHippocampalSubfields.sh bert right $inputDirectory $resultDirectory"

# Quantify the results, i.e., write a table for the left and right side with subfield volumes
doIt "kvlQuantifyHippocampalSubfieldSegmentations.sh $resultDirectory"


# Compare the values in these tables with those in reference tables, and count the number of failures
numberOfFailures=0;


# Loop over each side
sides="Left Right"
for side in $sides; do
  echo "Comparing side: $side"

  # Construct the table name and its reference for this brain side
  tableName="nonPartialVolumeStats"$side".txt"
  referenceTableName=""$atlasDirectory"/nonPartialVolumeStats"$side"ReferenceForBert.txt"

  # Read the contents of the tableName in an array, but skip the first line
  declare -a array
  array=( `tail -1 "$tableName"`)

  # Read the contents of the referenceTableName in an array, but skip the first line
  declare -a referenceArray
  referenceArray=( `tail -1 "$referenceTableName"`)

  # Loop over all the entries in the array, and compare each with its reference value
  startIndex=1  # Skip the first entry, which is "bert"
  let endIndex=${#array[*]}-1
  #for i in `seq $startIndex $endIndex`; do
  for i in `eval echo {$startIndex..$endIndex}`; do
    value=${array[ $i ]}
    referenceValue=${referenceArray[ $i ]}
    echo "  value: $value"
    echo "  referenceValue: $referenceValue"


    # Calculate the relative error
    float_scale=4  # Number of decimals behind the comma
    error=$(echo "scale=$float_scale; ( $value - $referenceValue) / $referenceValue" | bc -q 2>/dev/null )
    echo "  error: $error"

    # Take absolute value
    cond=$( echo "$error < 0" | bc -q 2>/dev/null)
    if [ $cond = 1 ]; then
      absError=$(echo "scale=$float_scale;  0 - $error"  | bc -q 2>/dev/null )
    else
      absError="$error"
    fi
    echo "  absError: $absError"

    # Test if it's small enough
    threshold="0.05" # 5% error is deemed acceptable
    cond=$( echo "$absError < $threshold" | bc -q 2>/dev/null)
    if [ $cond = 1 ]; then
      echo "  -> OK"
    else
      echo "  -> not OK"
      let numberOfFailures=$numberOfFailures+1
    fi
    echo "  --------"

  done # End loop over all volume entries

  echo "===================================="


done # End loop over left and right side


# Clean up
doIt "cd .."
doIt "rm -f $resultDirectory"



# Return the number of failures (can be checked by doing "echo $?" at the bash shell)
echo "numberOfFailures: $numberOfFailures"
exit $numberOfFailures






