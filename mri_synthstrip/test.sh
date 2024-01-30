#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# skull-stripped image
test_command mri_synthstrip -i in.mgz -o out.mgz
compare_vol out.mgz stripped.mgz

# GPU flags
test_command mri_synthstrip --gpu -g -i in.mgz --out out.mgz
compare_vol out.mgz stripped.mgz

# binary mask
test_command mri_synthstrip --image in.mgz -m out.mgz
compare_vol out.mgz mask.mgz

# distance transform without output validation
test_command mri_synthstrip -i in.mgz -d out.mgz
test_command mri_synthstrip -i in.mgz --sdt out.mgz

# default border value
test_command mri_synthstrip -b 1 -i in.mgz -m out.mgz
compare_vol out.mgz mask.mgz

# increased border
test_command mri_synthstrip --border 2 -i in.mgz -m out.mgz
compare_vol out.mgz border.mgz

# large border with SDT extension
test_command mri_synthstrip -b 8 -i in.mgz --mask out.mgz
compare_vol out.mgz large.mgz

# multiple frames, NIfTI format
test_command mri_synthstrip -i multi.in.nii.gz -m out.nii.gz
compare_vol out.nii.gz multi.mask.nii.gz
