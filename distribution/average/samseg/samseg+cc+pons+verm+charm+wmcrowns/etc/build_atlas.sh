#!/bin/bash

export FREESURFER_HOME=/autofs/cluster/koen/stefano/Block_to_Hemi/freesurfer/
source $FREESURFER_HOME/SetUpFreeSurfer.sh

atlas_dir=/autofs/cluster/koen/stefano/Block_to_Hemi/freesurfer/gems/bin/

data_dir=/autofs/cluster/fssubjects/atlases/fsv6.aseg_atlas/samseg-atlas/mni.22.12.31/

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

# Stiffness 0.1
# 20 iteration with edge collapse factor 1
# 20 iteration with edge collapse factor 1.005

python gems_train_atlas.py -n 3 -m 9 10 9 -s training-schedule.txt -l $data_dir/v6.981102_vc604.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.981112_vc623.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.981204_vc660.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990111_vc716.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990119_vc740.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990128_vc764.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990210_vc792.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990317_vc876.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990326_vc891.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990715_vc1131.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990729_vc1168.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.990903_vc1253.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991025_vc1379.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991102_vc1401.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991109_vc1420.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991110_vc1425.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991113_vc1439.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991113_vc1440.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991120_vc1456.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii $data_dir/v6.991122_vc1479.seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.nii -w /autofs/cluster/koen/stefano/Doug_Atlas/Test9109/

