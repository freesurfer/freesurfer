#!/usr/bin/env bash

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
    NNUNET_MODEL_DIR=$FREESURFER_HOME_FSPYTHON/models/nnUNetTrainer__nnUNetPlans__2d
fi

if [ ! -e "$NNUNET_SCRIPT_DIR" ]; then
    echo "WARNING: NNUNET_SCRIPT_DIR is not set"
    echo "Setting NNUNET_SCRIPT_DIR to be: $FREESURFER_HOME_FSPYTHON/python/packages/nnUNet_v2"
    NNUNET_SCRIPT_DIR=$FREESURFER_HOME_FSPYTHON/python/packages/nnUNet_v2
fi

# sanity check that the model files exist
ls "$NNUNET_MODEL_DIR/plans.json"
if [ $? -ne 0 ]; then
    echo "Unable to detect the required files in NNUNET_MODEL_DIR: $NNUNET_MODEL_DIR"
    echo "If this is your first time using this util, follow the instructions"
    echo "below to download, unpack, install, and configure the model files."
    echo ""
    echo "Download the model file from the FTP server. Copying the following link"
    echo "into your browser will begin the download: https://ftp.nmr.mgh.harvard.edu/pub/dist/lcnpublic/dist/dissection_photo_model/nnUNetTrainer__nnUNetPlans__2d.tar.gz"
    echo ""
    echo "Unpack the the downloaded tarball:"
    echo "  tar -xzf nnUNetTrainer__nnUNetPlans__2d.tar.gz"
    echo ""
    echo "Move the extracted directory to your models/ directory in your FS"
    echo "install. If FS was built with FSPYTHON_INSTALL_TREE=OFF, the models/"
    echo "dir will be located at \$FREESURFER_HOME/models"
    echo "If FSPYTHON_INSTALL_TREE=ON, then models/ will be located at"
    echo "\$FREESURFER_HOME/../fspython/python/models"
    echo ""
    echo "Once the model files are unpacked and installed into the models/dir,"
    echo "set NNUNET_MODEL_DIR to be the root of the directory you just moved"
    echo "into the models/ dir."
    echo ""
    echo "Then set NNUNET_SCRIPTS_DIR to be the path to the nnUNet_v2 dir under"
    echo "the packages/ directory."
    echo "If FSPYTHON_INSTALL_TREE=OFF was set during the build, this dir can"
    echo "be found under \$FREESURFER_HOME/python/packages."
    echo "If FSPYTHON_INSTALL_TREE=ON was set during the build, this dir can be"
    echo "under \$FREESURFER_HOME/../fspython/python/packages"
    exit 1
fi

fspython $NNUNET_SCRIPT_DIR/process_directory_no_upsample.py "$@"
