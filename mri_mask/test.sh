#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# untaring and taring the large gcam data is time consuming,
# so don't remove output before each test
FSTEST_NO_DATA_RESET=1 && init_testdata

# apply mask
test_command mri_mask nu.1.mgz brainmask.1.mgz mask.mgz
compare_vol mask.mgz mask.ref.mgz

# transform mask using LTA
test_command mri_mask -xform 2_to_1.lta nu.1.mgz brainmask.2.mgz mask.lta.mgz
compare_vol mask.lta.mgz mask.lta.ref.mgz

# transform mask using inverse LTA
test_command mri_mask -xform 1_to_2.lta nu.1.mgz brainmask.2.mgz mask.lta.inv.mgz
compare_vol mask.lta.inv.mgz mask.lta.inv.ref.mgz

# transform mask with FSLMAT
test_command mri_mask -xform 2_to_1.fslmat -lta_src nu.2.mgz -lta_dst nu.1.mgz nu.1.mgz brainmask.2.mgz mask.fsl.mgz
compare_vol mask.fsl.mgz mask.fsl.ref.mgz

# transform mask with GCAM
test_command mri_mask -xform 2_to_1.m3z nu.1.mgz brainmask.2.mgz mask.gcam.mgz
compare_vol mask.gcam.mgz mask.gcam.ref.mgz

# transfer edits
test_command mri_mask -transfer 255 -keep_mask_deletion_edits nu.2.mgz wm.2.mgz edits.mgz
compare_vol edits.mgz edits.ref.mgz

# transfer edits using identity.nofile
test_command mri_mask -xform identity.nofile -transfer 255 -keep_mask_deletion_edits nu.2.mgz wm.2.mgz edits.mgz
compare_vol edits.mgz edits.ref.mgz

# transfer edits using LTA
test_command mri_mask -xform 2_to_1.lta -transfer 255 -keep_mask_deletion_edits nu.1.mgz wm.2.mgz edits.lta.mgz
compare_vol edits.lta.mgz edits.lta.ref.mgz

# transfer edits using GCAM
test_command mri_mask -xform 2_to_1.m3z -transfer 255 -keep_mask_deletion_edits nu.1.mgz wm.2.mgz edits.gcam.mgz
compare_vol edits.gcam.mgz edits.gcam.ref.mgz

# transfer 255 only using inverse GCAM
test_command mri_mask -xform 1_to_2.m3z -transfer 255 nu.1.mgz wm.2.mgz edits.gcam.inv.mgz
compare_vol edits.gcam.inv.mgz edits.gcam.inv.ref.mgz
