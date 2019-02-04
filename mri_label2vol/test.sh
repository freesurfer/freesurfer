#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_label2vol --label label/lh.cortex.label --regheader mri/rawavg.mgz --temp mri/aparc+aseg.mgz --o l2v.nii
compare_vol l2v.nii l2v_ref.nii
