#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mris_curvature -seed 1234 -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated

compare_vol lh.inflated.H lh.inflated.ref.H
compare_vol lh.inflated.K lh.inflated.ref.K
