#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mris_sphere -seed 1234 rh.inflated rh.sphere
compare_surf rh.sphere rh.ref.sphere
