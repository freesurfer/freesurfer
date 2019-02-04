#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_seg_diff --debug --seg1 mri/aseg.auto.mgz --seg2 mri/aseg.mgz --diff aseg_diff_out.mgz
compare_vol aseg_diff_out.mgz aseg_diff_ref.mgz
