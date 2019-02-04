#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mri_coreg --mov template.nii.gz --targ orig.mgz --reg reg.lta --dof 12 --ftol .1 --linmintol .1
compare_file reg.lta source.lta -I#
