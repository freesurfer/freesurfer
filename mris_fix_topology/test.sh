#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 bert lh
compare_surf bert/surf/lh.orig bert/surf/lh.orig.ref