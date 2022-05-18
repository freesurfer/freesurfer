#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mri_watershed -T1 -brain_atlas ${FREESURFER_HOME}/average/RB_all_withskull_2016-05-10.vc700.gca talairach_with_skull.lta T1.mgz brainmask.mgz
compare_vol brainmask.mgz brainmask.ref${TESTDATA_SUFFIX}.mgz
