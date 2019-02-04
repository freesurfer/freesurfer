#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_ca_normalize \
    -mask brainmask.mgz \
    -c ctrl_pts.mgz \
    nu.mgz \
    ${FREESURFER_HOME}/average/RB_all_2016-05-10.vc700.gca \
    talairach.lta \
    norm.mgz

compare_vol norm.mgz norm.ref.mgz
compare_vol ctrl_pts.mgz ctrl_pts.ref.mgz
