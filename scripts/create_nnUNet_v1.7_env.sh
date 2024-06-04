#!/usr/bin/env bash

# This script will install the additional dependecies required for the 
# synthseg_tumor_nnUNet utility

### helper functions
# cleanup function to call on exit
cleanup(){
    echo "Cleaning up the tmp dir: $TMP_DIR"
    rm -rf $TMP_DIR
}
# so we can call 'realpath' on Mac 
realpath(){
    echo $(cd $(dirname $1); pwd)/$(basename $1)
}
# Show help text and exit with status 0|1
show_help(){
    echo "create_nnUNet_v1.7_env.sh"
    echo "      HELP"
    echo "USAGE:"
    echo "./create_nnUNet_v1.7_env.sh \ "
    echo "  -e <path to conda.sh> \ "
    echo "  -n <conda_env_name> \ "
    echo "  -m <path to Task002_StrokesLong>"
    echo "  -d <path to save nnUNet files>"
    echo ""
    echo "More detailed instructions can be found here under FREESURFER_HOME/docs/install_nnUNet_v1.7_instructions.md or on github: https://github.com/freesurfer/freesurfer/blob/dev/install_nnUNet_v1.7_instructions.md"
    exit $1
}

HELP=0
INSTALL_CUDA=0
CWD=$PWD

UNRECOGNIZED=()

# parse cli args
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            HELP=1
            shift
            ;;
        -c|--cuda)
            INSTALL_CUDA=1
            shift
            ;;
        -e|--conda-env)
            SOURCE_CONDA="$2"
            shift
            shift
            ;;
        -n|--name)
            CREATE_ENV_NAME="$2"
            shift
            shift
            ;;
        -m|--model-file)
            MODEL_FILE_PATH="$2"
            shift
            shift
            ;;
        -d|--nnUNet-model-destination)
            NNUNET_ENV_DIR="$2"
            shift
            shift
            ;;
        *)
            UNRECOGNIZED+=("$1")
            shift
            ;;
    esac
done

### validate args
# help check
if [[ $HELP -eq 1 ]]; then
    show_help 0
fi

# check for unrecognized args
if [[ ${#UNRECOGNIZED} -ne 0 ]]; then
    echo "Unrecognized command line args:"
    echo "${UNRECOGNIZED[@]}"
    show_help 1
fi

# test that the conda.sh script exists
if [ ! -f $SOURCE_CONDA ]; then
    echo "Cannot access the 'conda.sh' script at $SOURCE_CONDA, exiting..."
    show_help
fi

# test that the model file exists
if [ ! -d $MODEL_FILE_PATH ]; then
    echo "The directory containing the model files does not exist: $MODEL_FILE_PATH"
    show_help
fi

# test that the path to store the nnUNet file tree is set
if [ -z $NNUNET_ENV_DIR ]; then
    echo "The directory to save the nnUNet model files is not set. exiting..."
    show_help
fi

### build paths for the nnUNet directory structure
#set the path to the root of the nnUNet file tree
MODEL_DESTINATION=$NNUNET_ENV_DIR/nnUNet_v1.7

# directory structure required by nnUNet
# need to create these then call realpath, or the Mac hack breaks
nnUNet_raw_data_base=$MODEL_DESTINATION/nnUNet_raw_data_base
nnUNet_raw_data=$nnUNet_raw_data_base/nnUNet_raw_data
nnUNet_preprocessed=$MODEL_DESTINATION/nnUNet_preprocessed
RESULTS_FOLDER=$MODEL_DESTINATION/nnUNet_trained_models
TRAINED_MODEL_DESTINATION=$RESULTS_FOLDER/nnUNet/3d_fullres

# create the file hierarchy needed by nnUNet
mkdir -p $MODEL_DESTINATION $nnUNet_raw_data $nnUNet_preprocessed $RESULTS_FOLDER $TRAINED_MODEL_DESTINATION

nnUNet_raw_data_base=`realpath $MODEL_DESTINATION/nnUNet_raw_data_base`
nnUNet_preprocessed=`realpath $MODEL_DESTINATION/nnUNet_preprocessed`
RESULTS_FOLDER=`realpath $MODEL_DESTINATION/nnUNet_trained_models`
TRAINED_MODEL_DESTINATION=`realpath $RESULTS_FOLDER/nnUNet/3d_fullres`

if [[ $? -ne 0 ]]; then
    echo "An error occurred while attempting to create the nnUNet directories"
    exit 1
fi

# copy over the model file
cp -r $MODEL_FILE_PATH $TRAINED_MODEL_DESTINATION
if [[ $? -ne 0 ]]; then
    echo "An error occurred while attempting to copy the trained model files."
    exit 1
fi

### tmp dir setup
# create the tmp dir
TMP_DIR=$(mktemp -d)

# set trap to call cleanup on exit
trap cleanup EXIT

### source conda and build the env
# source the conda.sh script
source $SOURCE_CONDA

# cretate the conda env
CONDA_CREATE_CMD="conda create -n $CREATE_ENV_NAME python=3.8.13 pytorch=2.1.2 torchvision=0.16.2 torchaudio=2.1.2 -c pytorch -y"

# append the cuda dependencies if --cuda passed
if [[ $INSTALL_CUDA -eq 1 ]]; then
    CONDA_CREATE_CMD="$CONDA_CREATE_CMD pytorch-cuda=11.8 -c nvidia"
fi

echo "Creating conda env: $CONDA_CREATE_CMD"
exit 0
$CONDA_CREATE_CMD
if [[ $? -ne 0 ]]; then 
    echo "An error occurred while trying to create the conda env."
    exit 1
fi

# activate the conda env
conda activate $CREATE_ENV_NAME

if [[ $? -ne 0 ]]; then
    echo "An error occurred attempting to activate the conda env: $CONDA_ENV_NAME"
    echo "Delete that environment before attempting to run this script again"
    exit 1
fi

### nnUNet download and install
# git clone the nnUNet source, check out the proper commit, pip install it
GIT_CLONE_CMD="git clone https://github.com/MIC-DKFZ/nnUNet.git $TMP_DIR/nnUNet"
GIT_COMMIT_HASH="7f1e273fa1021dd2ff00df2ada781ee3133096ef"
GIT_CHECKOUT_CMD="git checkout $GIT_COMMIT_HASH"

echo "Cloning nnUNet source code"
$GIT_CLONE_CMD
if [[ $? -ne 0 ]]; then 
    echo "An error occurred cloning the nnUNet source code."
    exit 1
fi

# cd into the git repo to checkout the commit
echo "Checking out the proper commit"
cd $TMP_DIR/nnUNet
$GIT_CHECKOUT_CMD

if [[ $? -ne 0 ]]; then
    echo "An error occurred checking out the specific nnUNet commit"
    exit 1
fi

# pip install nnUNet source
pip install .
if [[ $? -ne 0 ]]; then
    echo "An error occurred installing the nnUNet python module"
    exit 1
fi

# cd back to the old PWD
cd $CWD

PATH_FILE="nnUNet_v1.7_path.sh"
# add bash shebang to PATH_FILE
echo "#!/usr/bin/env bash" >> $PATH_FILE

echo ""
echo "Set the following variables in your env so nnUNet can locate your models and data"
echo "export nnUNet_raw_data_base=$nnUNet_raw_data_base" | tee -a $PATH_FILE
echo "export nnUNet_preprocessed=$nnUNet_preprocessed" | tee -a $PATH_FILE
echo "export RESULTS_FOLDER=$RESULTS_FOLDER" | tee -a $PATH_FILE

echo "Installation complete"
echo "Activate the conda envrionment named: $CREATE_ENV_NAME, and source $PATH_FILE to run the nnUNet model."
echo "You can save a copy of $PATH_FILE in $CONDA_PREFIX/etc/conda/activate.d, for this to happen automatically when you activate the environment."
