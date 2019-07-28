#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# untaring and taring the large gcam data is time consuming,
# so don't remove output before each test
FSTEST_NO_DATA_RESET=1 && init_testdata

# copy LTA
test_command mri_concatenate_gcam im1_to_im2.lta copy.lta
compare_lta copy.lta copy.ref.lta

# concatenate LTAs
test_command mri_concatenate_gcam im3_to_im1.lta im1_to_im2.lta concat.lta
compare_lta concat.lta concat.ref.lta

# make sure we have mri_convert in the build tree before continuing
mri_convert=$(find_path $FSTEST_CWD mri_convert/mri_convert)

# change the GCAM source image space
test_command mri_concatenate_gcam -s im1.mgz im2_to_im3.m3z gcam.m3z
test_command $mri_convert -at gcam.m3z im2_in_space1.mgz src.mgz
compare_vol src.mgz src.ref.mgz

# change the GCAM target image space
test_command mri_concatenate_gcam -t im1.mgz im2_to_im3.m3z gcam.m3z
test_command $mri_convert -at gcam.m3z im2.mgz dest.mgz
compare_vol dest.mgz dest.ref.mgz

# change the GCAM target image space to a larger space
test_command mri_concatenate_gcam --change-target im2.mgz im2_to_im1.m3z gcam.m3z
compare_vol gcam.m3z larger.ref.m3z

# concatenate LTA with GCAM
test_command mri_concatenate_gcam im1_to_im2.lta im2_to_im3.m3z gcam.m3z
test_command $mri_convert -at gcam.m3z im1.mgz lta_gcam.mgz
compare_vol lta_gcam.mgz lta_gcam.ref.mgz

# concatenate GCAM with LTA
test_command mri_concatenate_gcam im2_to_im3.m3z im3_to_im1.lta gcam.m3z
test_command $mri_convert -at gcam.m3z im2.mgz gcam_lta.mgz
compare_vol gcam_lta.mgz gcam_lta.ref.mgz

# concatenate LTA with GCAM with LTA
test_command mri_concatenate_gcam im1_to_im2.lta im2_to_im3.m3z im3_to_im1.lta gcam.m3z
test_command $mri_convert -at gcam.m3z im1.mgz lta_gcam_lta.mgz
compare_vol lta_gcam_lta.mgz lta_gcam_lta.ref.mgz

# concatenate two GCAMs
test_command mri_concatenate_gcam im2_to_im3.m3z im3_to_im2.m3z gcam.m3z
test_command $mri_convert -at gcam.m3z im2.mgz gcam_gcam.mgz
compare_vol gcam_gcam.mgz gcam_gcam.ref.mgz

# downsample GCAM - leave some capacity for change
test_command mri_concatenate_gcam -d im2_to_im3.m3z gcam.m3z
test_command $mri_convert -at gcam.m3z im2.mgz down.mgz
compare_vol down.mgz down.ref.mgz --thresh 0.5
