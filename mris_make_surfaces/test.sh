#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

test_command mris_make_surfaces -aseg aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs subject lh
compare_surf subject/surf/lh.white.preaparc subject/surf/lh.white.preaparc.REF

test_command mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg aseg.presurf -mgz -T1 brain.finalsurfs subject lh
compare_surf subject/surf/lh.white subject/surf/lh.white.REF
compare_surf subject/surf/lh.pial subject/surf/lh.pial.REF
