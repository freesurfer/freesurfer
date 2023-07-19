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

if [ "$host_os" == "macos12" ]; then
   TESTDATA_SUFFIX=".clang13"
fi

if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "macos10" ]] || [[ "$host_os" == "macos12" ]] ; then
   compare_vol talairach.m3z talairach.ref${TESTDATA_SUFFIX}.m3z
else
   compare_vol talairach.m3z talairach.ref.m3z
fi

