#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_ca_register \
    -nobigventricles \
    -T talairach.lta \
    -align-after \
    -levels 3 \
    -n 2 \
    -tol 1.0 \
    -mask brainmask.mgz \
    norm.mgz \
    ${FREESURFER_HOME}/average/RB_all_2016-05-10.vc700.gca \
    talairach.m3z

compare_vol talairach.m3z talairach.ref.m3z
