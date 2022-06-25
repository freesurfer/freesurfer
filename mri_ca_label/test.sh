#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz talairach.m3z \
    ${FREESURFER_HOME}/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.out.mgz

if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "macos10" ]]; then
   compare_vol aseg.auto_noCCseg.out.mgz aseg.auto_noCCseg${TESTDATA_SUFFIX}.mgz
else
   compare_vol aseg.auto_noCCseg.out.mgz aseg.auto_noCCseg.mgz
fi

