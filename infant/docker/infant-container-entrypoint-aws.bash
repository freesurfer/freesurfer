#!/bin/bash

echo "==============================================================="
echo "ENVIRONMENT VARIABLES"
echo ""
echo "AWS_BATCH_JOB_ID:                  $AWS_BATCH_JOB_ID"
echo "AWS_BATCH_JQ_NAME:                 $AWS_BATCH_JQ_NAME"
echo "AWS_BATCH_CE_NAME:                 $AWS_BATCH_CE_NAME"
echo "---------------------------------------------------------------"
echo "FREESURFER_HOME:                   $FREESURFER_HOME"
echo "FS_INFANT_MODEL:                   $FS_INFANT_MODEL"
echo "SUBJECTS_DIR:                      $SUBJECTS_DIR"
echo "FS_SUB_NAME:                       $FS_SUB_NAME"
echo "---------------------------------------------------------------"
echo "FS_NIFTI_INPUT_S3_FILEPATH:        $FS_NIFTI_INPUT_S3_FILEPATH"
echo "FS_NIFTI_INPUT_LOCAL_FILEPATH:     $FS_NIFTI_INPUT_LOCAL_FILEPATH"
echo "FS_OUTPUT_S3_FILEPATH:             $FS_OUTPUT_S3_FILEPATH"
echo "==============================================================="

# infant pipeline input must be ${SUBJECTS_DIR}:${FS_SUB_NAME}/mprage.nii.gz
if [ -n "${SUBJECTS_DIR}" ] && [ -n "${FS_SUB_NAME}" ]; then
  echo "---------------------------------------------------------------"
  echo "SUBJECTS_DIR and FS_SUB_NAME detected. Attempting to make dir"
  echo "mkdir -p ${SUBJECTS_DIR}/${FS_SUB_NAME}"
  mkdir -p ${SUBJECTS_DIR}/${FS_SUB_NAME}
  echo "---------------------------------------------------------------"
fi

if [ -n "${FS_NIFTI_INPUT_S3_FILEPATH}" ] && [ -n "${FS_NIFTI_INPUT_LOCAL_FILEPATH}" ]; then
  echo "---------------------------------------------------------------"
  echo "FS_NIFTI_INPUT_S3_FILEPATH and FS_NIFTI_INPUT_LOCAL_FILEPATH detected. Attempting to copy file locally"
  echo "aws s3 cp ${FS_NIFTI_INPUT_S3_FILEPATH} ${FS_NIFTI_INPUT_LOCAL_FILEPATH}"
  aws s3 cp $FS_NIFTI_INPUT_S3_FILEPATH $FS_NIFTI_INPUT_LOCAL_FILEPATH
  echo "---------------------------------------------------------------"
fi

# Symlink the volume_mounted model files to where FreeSurfer expects them
# -----------------------------------------------------------------------
mkdir -p $FREESURFER_HOME/average/sscnn_skullstripping
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/cor_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/cor_sscnn.h5
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/ax_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/ax_sscnn.h5
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/sag_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/sag_sscnn.h5 

eval "$@"

if [ -n "${FS_OUTPUT_S3_FILEPATH}" ]; then
  echo "---------------------------------------------------------------"
  echo "FS_OUTPUT_S3_FILEPATH detected. Attempting to copy the subjects_dir to s3:"
  echo "aws s3 cp --recursive ${SUBJECTS_DIR}/${FS_SUB_NAME} ${FS_OUTPUT_S3_FILEPATH}"
  aws s3 cp --recursive ${SUBJECTS_DIR}/${FS_SUB_NAME}  ${FS_OUTPUT_S3_FILEPATH}
  echo "---------------------------------------------------------------"
fi
