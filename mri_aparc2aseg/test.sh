#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_aparc2aseg --s bert
compare_vol bert/mri/aparc+aseg.mgz bert/mri/aparc+aseg.ref.mgz

test_command mri_aparc2aseg --s bert --labelwm --hypo-as-wm --rip-unknown --o ${SUBJECTS_DIR}/bert/mri/wmparc.mgz
compare_vol bert/mri/wmparc.mgz bert/mri/wmparc.ref.mgz
