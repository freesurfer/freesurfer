#! /bin/bash

# check to see that MCR is installed
MCRstatus=eval checkMCR.sh

# if MCR is not found, will print instructions on how to
# install MCR
if [[ ! -z $MCRstatus ]]; then
    echo $MCRstatus
    exit 1
fi


# parser args
POSITIONAL_ARGS=()
# show help flag
SHOW_HELP=0
while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--subject)
      SUBJECT="$2"
      shift
      shift
      ;;
    -d|--dti_dir)
      DTI_DIR="$2"
      shift
      shift
      ;;
    -o|--out_dir)
      OUT_DIR="$2"
      shift
      shift
      ;;
    --hrdti)
      HRDTI="$2"
      shift
      shift
      ;;
    --use_unregistered)
      USE_UNREGISTERED="$2"
      shift
      shift
      ;;
    -h|--help)
      SHOW_HELP=1
      shift
      ;;
  esac
done

# if --h|--help passed, print help text and exit
if [ $SHOW_HELP == 1 ]; then
    echo " "
    echo "USAGE:"
    echo "segmentThalamicNuclei.sh -s [SUBJECT_ID]"
    echo "segmentThalamicNuclei.sh -s [SUBJECT_ID] -o [OUT_DIR] -d [DTI_DIR] --hrdti [0|1] --use_unregistered [0|1]"
    echo " "
    echo "More information can be found at: https://surfer.nmr.mgh.harvard.edu/fswiki/ThalamicNucleiDTI"
    exit
fi

# show command line args and important env vars
set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters
echo "CLI ARGS:"
echo "SUBJECT           = ${SUBJECT}"
echo "DTI_DIR           = ${DTI_DIR}"
echo "OUT_DIR           = ${OUT_DIR}"
echo "HRDTI             = ${HRDTI}"
echo "USE_UNREGISTERED  = ${USE_UNREGISTERED}"
echo " "
echo "ENV VARS:"
echo "FREESURFER_HOME   = ${FREESURFER_HOME}"
echo "SUBJECTS_DIR      = ${SUBJECTS_DIR}"
echo "THALSEG_HOME      = ${THALSEG_HOME}"

# ITK THREADS
if [[ -z $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS ]]; then
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
elif [ $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS == 0 ]; then
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
fi

# is_running_file/log file
IS_RUNNING_NAME="segmentThalamicNuclei_DTI_ISRUNNING.txt"
LOG_NAME="segmentThalamicNuclei_DTI.log"
# if outdir is set, put the is_running file there
if [[ ! -z $OUT_DIR ]]; then
    IS_RUNNING_FILE=$OUT_DIR/$IS_RUNNING_NAME
    LF=$OUT_DIR/$LOG_NAME
else
    IS_RUNNING_FILE=$SUBJECTS_DIR/$SUBJECT/$IS_RUNNING_NAME
    LF=$SUBJECTS_DIR/$SUBJECT/$LOG_NAME
fi

# check if the is_running_file exists
if [[ -e $IS_RUNNING_FILE ]]; then
    echo "$SUBJECTS_DIR/$SUBJECT is actively running, allow it to finish, or run:"
    echo "rm $SUBJECTS_DIR/$SUBJECT/$IS_RUNNING_NAME"
    exit
fi

# Set the ref to MCR and build the command to run the actual exe wrapper
FS_MCR=$FREESURFER_HOME/MCRv97

# building the command, adding in --flags is annoying, that's why they compiled with positional args
cmd="run_SegmentThalamicNuclei_DTI.sh $FS_MCR -s $SUBJECT"

if [[ ! -z $OUT_DIR ]]; then
    cmd="$cmd -o $OUT_DIR"
fi

if [[ ! -z $DTI_DIR ]]; then
    cmd="$cmd -d $DTI_DIR"
fi

#if [ $HRDTI == 1 ]; then
    # concat flag for --hrdti
#    echo ""
#fi

#if [ $USE_UNREGISTERED == 1 ]; then
#    # concat flag for --use_unregistered
#    echo ""
#fi

echo "COMMAND:"
echo "$cmd" | tee -a $LF

touch $IS_RUNNING_FILE

$cmd | tee -a $LF

rm $IS_RUNNING_FILE

echo "If used for publications, please site:"
echo "Accurate Bayesian segmentation of thalamic nuclei using diffusion MRI and an improved histological atlas. "
echo "Tregidgo HFJ, Soskic S, Althonayan J, Maffei C, Van Leemput K, Golland P, Insausti R, Lerma-Usabiaga G, "
echo "Caballero-Gaudes C, Paz-Alonso PM, Yendiki A, Alexander DC, Bocchetta M, Rohrer JD, Iglesias JE., NeuroImage, 274, 120129,(2023)"