#!/bin/bash 

# check for required env vars
if [ ! -e "$FREESURFER_HOME" ]; then
    echo "error: FREESURFER_HOME has not been properly set"
    exit 1
fi

if [ ! -e "$NNUNET_MODEL_DIR" ]; then
    echo "WARNING: NNUNET_MODEL_DIR is not set"
    echo "This should be the path to the directory containing the model"
    echo "dependency 'nnUNetTrainer__nnUNetPlans__2d'"
    echo "Setting NNUNET_MODEL_DIR to be: $FREESURFER_HOME_FSPYTHON/models/nnUNetTrainer__nnUNetPlans__2d"
    export NNUNET_MODEL_DIR=$FREESURFER_HOME_FSPYTHON/models/nnUNetTrainer__nnUNetPlans__2d
fi

if [ ! -e "$NNUNET_SCRIPT_DIR" ]; then
    echo "WARNING: NNUNET_SCRIPT_DIR is not set"
    echo "Setting NNUNET_SCRIPT_DIR to be: $FREESURFER_HOME_FSPYTHON/python/packages/nnUNet_v2"
    export NNUNET_SCRIPT_DIR=$FREESURFER_HOME_FSPYTHON/python/packages/nnUNet_v2
fi

# sanity check that the model files exist
ls "$NNUNET_MODEL_DIR/plans.json"
if [ $? -ne 0 ]; then
    echo "Unable to detect the required files in NNUNET_MODEL_DIR: $NNUNET_MODEL_DIR"
    echo "If this is your first time using this util, follow the instructions"
    echo "below to download, unpack, install, and configure the model files."
    echo ""
    echo "  1. Download the model file from the FTP server. Copying the following link"
    echo "     into your browser will begin the download: https://ftp.nmr.mgh.harvard.edu/pub/dist/lcnpublic/dist/dissection_photo_model/nnUNetTrainer__nnUNetPlans__2d.tar.gz"
    echo ""
    echo "  2. Move the tarball to the FreeSurfer modles/ dir and cd into it:"
    echo "     mv nnUNetTrainer__nnUNetPlans__2d.tar.gz \$FREESURFER_HOME/models"
    echo "     cd \$FREESURFER_HOME/models"
    echo ""
    echo "  3. Unpack the tarball:"
    echo "     tar -xzf nnUNetTrainer__nnUNetPlans__2d.tar.gz"
    echo ""
    echo "  4. Configure the nnUNet environment variables:"
    echo "     export NNUNET_MODEL_DIR=\$FREESURFER_HOME/models/nnUNetTrainer__nnUNetPlans__2d"
    echo "     export NNUNET_SCRIPTS_DIR=\$FREESURFER_HOME/python/packages/nnUNet_v2"
    echo ""
    exit 1
fi

dissection_photo "$@"
