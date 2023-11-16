#!/bin/bash

# USAGE: ./quantifyData.sh -f <data_file.txt> [-f <data_file_2.txt>] -o out_dir [-s SUBJECTS_DIR]
# Will concatenate stats for all subjects found in SUBJECTS_DIR and output them
# to out_dir, a file will be generated for each data file specified with the -f
# flag. Each output file will have the naming convention <data_file>_concat.txt
# SUBJECTS_DIR env var will be respected or can be overwritten by specifying -s
# If out_dir is not specified, results files will output to SUBJECTS_DIR

# Define a function that exits if something goes wrong.
function doIt {

  command="$1"

  eval "$command"

  if [ $? != 0 ]; then
    echo "failed to do $command"
    exit -1
  fi
}

# echo the original command line
echo "cmd: $0 $*"

# help flag
HELP=0

if [[ $# -le 2 ]]; then
    HELP=1
fi

# parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--subject)
            SUBJECTS_DIR=$2
            shift
            shift
            ;;
        -o|--output)
            # add a '/' to help with path concat later
            OUTPATH="$2/"
            shift
            shift
            ;;
        -f|--file)
            FILES+=("$2")
            shift
            shift
            ;;
        -h|--help)
            HELP=1
            shift
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# check if help flag passed, print usage and exit
if [[ $HELP -eq 1 ]]; then
    echo ""
    echo "USAGE: $0 -f <data_file.txt> [-f <data_file_2.txt>] [-o out_dir] [-s SUBJECTS_DIR]"
    echo ""
    echo "Will concatenate stats for all subjects found in SUBJECTS_DIR and output them"
    echo "to out_dir, a file will be generated for each data file specified with the -f"
    echo "flag. Each output file will have the naming convention <data_file>_concat.txt"
    echo "SUBJECTS_DIR env var will be respected or can be overwritten by specifying -s"
    echo "If out_dir is not specified, results files will output to the SUBJECTS_DIR"
    exit
fi

# echo parsed args
echo "Gathering results from: "
echo "$SUBJECTS_DIR"
echo "Looking at files:"
echo "${FILES[@]}"
echo "Outputting results to:"
if [[ -z $OUTPATH ]]; then
    echo "$SUBJECTS_DIR"
else
    echo "$OUTPATH"
fi

# echo any unused args if passed
if [[ ${#POSITIONAL[@]} -ne 0 ]]; then
    echo "Unrecognized args:"
    echo "${POSITIONAL[@]}"
fi

# go to subjects dir
doIt "cd $SUBJECTS_DIR"

# get list of subjects to parse
SUB_NAMES=(`ls -d */`)
NUM_SUBS=${#SUB_NAMES[@]}

# concatenate stats for all files specified
for stat_file in ${FILES[@]}; do
    # set outfile name for related stat_file
    results_file="$OUTPATH${stat_file%'.txt'}_concat.txt"
    echo ""
    echo "Concatenating stats for file: $stat_file"
    echo "Will write results to: $results_file"


    WRITE_HEADER='yes'
    for subject in ${SUB_NAMES[@]}; do
        # subject name
        subjectName=`echo "${subject//\/}"`
        echo "Working on subject: $subjectName"

        # volume files tp concat
        VolFile="$subjectName/mri/$stat_file"

        # check that the VolFiles exist
        if [ -f $VolFile ]; then
        # If first subjet, write file header
            if [ $WRITE_HEADER == 'yes' ]; then
                WRITE_HEADER='no'

                header_string="Subject"
                # add left labels to header
                while read line; do
                    arr=(`echo ${line}`)
                    header_string="$header_string left_${arr[0]}"
                done < $VolFile

                #write header to file
                echo $header_string > $results_file
            fi

            # Get data from Vol files and append to outfile
            data_string="$subjectName"
            # add left hemi data
            while read line; do
                arr=(`echo ${line}`)
                data_string="$data_string ${arr[1]}"
            done < $VolFile
        
            # write data to file
            echo $data_string >> $results_file
        else
            echo "Could not find $stat_file, skipping subject."
            echo ""
        fi
    done
done