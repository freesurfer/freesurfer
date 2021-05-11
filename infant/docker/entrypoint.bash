#!/bin/bash

# Symlink the volume_mounted model files to where FreeSurfer expects them
mkdir -p $FREESURFER_HOME/average/sscnn_skullstripping
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/cor_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/cor_sscnn.h5
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/ax_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/ax_sscnn.h5
ln -s $FS_INFANT_MODEL/sscnn_skullstrip/sag_sscnn.h5 $FREESURFER_HOME/average/sscnn_skullstripping/sag_sscnn.h5 

eval "$@"
