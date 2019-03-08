#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

# modified (shortened) usage in recon-all
test_command mris_register -curv -rusage rusage.mris_register.lh.dat lh.sphere lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif lh.sphere.reg
compare_surf lh.sphere.reg ref_lh.sphere.reg
