#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# avoid untaring before each test to save time
FSTEST_NO_DATA_RESET=1 && init_testdata

# LPS-to-LPS displacement field to M3Z
test_command mri_warp_convert --inlps lps.nii.gz --insrcgeom geom.mgz --outm3z out.m3z
compare_vol out.m3z ref.m3z --thresh 0.001

# RAS-to-RAS displacement field to M3Z
test_command mri_warp_convert --inras ras.nii.gz --insrcgeom geom.mgz --outm3z out.m3z
compare_vol out.m3z ref.m3z --thresh 0.001

# M3Z to LPS-to-LPS deplacement field
test_command mri_warp_convert --inm3z ref.m3z --outlps out.nii.gz
compare_vol out.nii.gz lps.nii.gz --thresh 0.001

# M3Z to RAS-to-RAS deplacement field
test_command mri_warp_convert --inm3z ref.m3z --outras out.nii.gz
compare_vol out.nii.gz ras.nii.gz --thresh 0.001

# LPS-to-LPS to RAS-to-RAS
test_command mri_warp_convert --inlps lps.nii.gz --insrcgeom geom.mgz --outras out.nii.gz
compare_vol out.nii.gz ras.nii.gz --thresh 0.001

# RAS-to-RAS to LPS-to-LPS
test_command mri_warp_convert --inras lps.nii.gz --insrcgeom geom.mgz --outras out.nii.gz
compare_vol out.nii.gz lps.nii.gz --thresh 0.001

# downsample M3Z
test_command mri_warp_convert --inlps lps.nii.gz --insrcgeom geom.mgz --outm3z out.m3z --downsample
compare_vol out.m3z half.ref.m3z --thresh 0.001

# downsample M3Z
test_command mri_warp_convert --inras ras.nii.gz -g geom.mgz --outm3z out.m3z -d
compare_vol out.m3z half.ref.m3z --thresh 0.001

# convert downsampled M3Z
test_command mri_warp_convert --inm3z half.ref.m3z --outlps out.nii.gz
compare_vol out.nii.gz half.ref.nii.gz --thresh 0.001

# LPS to source-voxel shift
test_command mri_warp_convert --inlps lps.nii.gz -g geom.mgz --outvox out.mgz
compare_vol out.mgz vox.mgz --thresh 0.001

# source-voxel shift to LPS
test_command mri_warp_convert --invox vox.mgz -g geom.mgz --outlps out.nii.gz
compare_vol out.nii.gz lps.nii.gz --thresh 0.001

# convert m3z to 4D mgz warp
test_command mri_warp_convert --inm3z ref.m3z --outmgzwarp out.m3z.mgz --outwarpformat abs-crs
compare_vol out.m3z.mgz ref.m3z.mgz
