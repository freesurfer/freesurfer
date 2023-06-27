#!/bin/bash

export FREESURFER_HOME=/usr/local/freesurfer/dev/
source $FREESURFER_HOME/SetUpFreeSurfer.sh

mesh_collection="/autofs/cluster/koen/stefano/Doug_Atlas/Test9109/epoch_0001_out_dir/CurrentMeshCollection20.gz"
mesh_level_1="/autofs/cluster/koen/stefano/Doug_Atlas/Test9109/epoch_0002_out_dir/CurrentMeshCollection5.gz"
mesh_affine="/autofs/cluster/koen/stefano/Doug_Atlas/Test9109/epoch_0002_out_dir/CurrentMeshCollection20.gz"
compression_lut="/autofs/cluster/koen/stefano/Doug_Atlas/compression_lookup_table_correct_names.txt" # LUT should automatically import names from another lut
template="/autofs/cluster/koen/stefano/Doug_Atlas/template_mask.nii.gz"  # TODO: this should be created automatically!
atlas_out_dir="/autofs/cluster/koen/stefano/Doug_Atlas/SAMSEG_atlas_folder/"
gmm_file_name="/autofs/cluster/koen/stefano/Doug_Atlas/shared_GMM_parameter_file.txt"
fs_lookup_table="/autofs/cluster/fssubjects/atlases/fsv6.aseg_atlas/samseg-atlas/mni.22.12.31/FreeSurferColorLUT3.txt"


fspython /autofs/cluster/koen/stefano/Block_to_Hemi/freesurfer/python/gems/prepareAtlasDirectory.py --mesh $mesh_collection --compression_lut $compression_lut --template $template --template $template --atlas $atlas_out_dir --shared_gmm_params $gmm_file_name --mesh_level_1 $mesh_level_1 --mesh_affine $mesh_affine --fs_lookup_table $fs_lookup_table -ac Unknown -ac GlobalWM -ac Veins GlobalGM Putamen Pallidum GlobalCSF -ac Tissues Spine Spinal-Cord
