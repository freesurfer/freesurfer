#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

# make surfaces

# test_command mris_make_surfaces -aseg aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs subject lh
# compare_surf subject/surf/lh.white.preaparc subject/surf/lh.white.preaparc.REF

test_command mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg aseg.presurf -mgz -T1 brain.finalsurfs subject lh
compare_surf subject/surf/lh.white subject/surf/lh.white.REF
compare_surf subject/surf/lh.pial subject/surf/lh.pial.REF

# place surfaces

export SUBJECTS_DIR=$SUBJECTS_DIR/subject
test_command "cd $SUBJECTS_DIR/mri && mris_place_surface --adgws-in ../surf/autodet.gw.stats.lh.dat --wm wm.mgz --threads 1 --invol brain.finalsurfs.mgz --lh --i ../surf/lh.orig --o ../surf/lh.white.preaparc --white --seg aseg.presurf.mgz --nsmooth 5"

export SUBJECTS_DIR=$SUBJECTS_DIR
# compare_surf surf/lh.white.preaparc surf/lh.white.preaparc.REF_PLACE_SURFACES
(cd $SUBJECTS_DIR && ../../../mris_diff/mris_diff surf/lh.white.preaparc surf/lh.white.preaparc.REF_PLACE_SURFACES)

