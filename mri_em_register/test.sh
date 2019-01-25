#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mri_em_register -uns 3 -mask brainmask.mgz nu.mgz ${FREESURFER_HOME}/average/RB_all_2016-05-10.vc700.gca talairach.lta
compare_lta talairach.lta talairach.ref.lta
