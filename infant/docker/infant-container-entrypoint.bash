#!/bin/bash

# TODO: remove `FS_INFANT_MODEL` (no longer needed?)
echo "==============================================================="
echo "ENVIRONMENT VARIABLES"
echo ""
echo "FREESURFER_HOME:                   $FREESURFER_HOME"
echo "FS_INFANT_MODEL:                   $FS_INFANT_MODEL"
echo "SUBJECTS_DIR:                      $SUBJECTS_DIR"
echo "FS_SUB_NAME:                       $FS_SUB_NAME"
echo "SSCNN_MODEL_DIR:                   $SSCNN_MODEL_DIR"
echo "==============================================================="

# Symlink the volume_mounted model files to where FreeSurfer expects them
# PW 2021/11/18 No longer needed, since `SSCNN_MODEL_DIR` can now be used in `sscnn_skullstrip`
# -----------------------------------------------------------------------
#mkdir -p $FREESURFER_HOME/average/sscnn_skullstripping
#ln -s $FS_INFANT_MODEL/sscnn_skullstrip/cor_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/cor_sscnn.h5
#ln -s $FS_INFANT_MODEL/sscnn_skullstrip/ax_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/ax_sscnn.h5
#ln -s $FS_INFANT_MODEL/sscnn_skullstrip/sag_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/sag_sscnn.h5 

eval "$@"
