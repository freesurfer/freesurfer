#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_warp_convert --inm3z talairach.m3z --outitk out.nii.gz
compare_vol out.nii.gz ref.nii.gz
