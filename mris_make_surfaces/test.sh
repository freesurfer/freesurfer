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

# cmd 1
#    input = lh_DECIMATE_AREA_5.orig
#    output = ./surf/lh_DECIMATE_AREA_5.white.preaparc  

cmd_1_banner="========== COMMAND 1 for mris_place_surfaces =========="
cmd_1="mris_place_surface --adgws-in ../surf/autodet.gw.stats.lh_DECIMATE_AREA_5.dat --wm wm.mgz --threads 1 --invol brain.finalsurfs_DECIMATE_AREA_5.mgz --lh --i ../surf/lh_DECIMATE_AREA_5.orig --o ../surf/lh_DECIMATE_AREA_5.white.preaparc --white --seg aseg.presurf_DECIMATE_AREA_5.mgz --nsmooth 5"

# cmd 2
#    input  = ./surf/lh.white_DECIMATE_AREA_5.preaparc (output from cmd 1) 
#    output = ./surf/lh_DECIMATE_AREA_5.white

cmd_2_banner="========== COMMAND 2 for mris_place_surfaces =========="
cmd_2="mris_place_surface --adgws-in ../surf/autodet.gw.stats.lh_DECIMATE_AREA_5.dat --seg aseg.presurf_DECIMATE_AREA_5.mgz --threads 1 --wm wm_DECIMATE_AREA_5.mgz --invol brain.finalsurfs_DECIMATE_AREA_5.mgz --lh --i ../surf/lh_DECIMATE_AREA_5.white.preaparc --o ../surf/lh_DECIMATE_AREA_5.white --white --nsmooth 0 --rip-label ../label/lh_DECIMATE_AREA_5.cortex.label --rip-bg --rip-surf ../surf/lh_DECIMATE_AREA_5.white.preaparc --aparc ../label/lh_DECIMATE_AREA_5.aparc.annot"

# cmd 3
#    input = ./surf/lhi_DECIMATE_AREA_5.white (output from cmd 2)
#    output = ./surf/lh_DECIMATE_AREA_5.pial.T1

cmd_3_banner="========== COMMAND 3 for mris_place_surfaces =========="
cmd_3="mris_place_surface --adgws-in ../surf/autodet.gw.stats.lh_DECIMATE_AREA_5.dat --seg aseg.presurf_DECIMATE_AREA_5.mgz --threads 1 --wm wm_DECIMATE_AREA_5.mgz --invol brain.finalsurfs_DECIMATE_AREA_5.mgz --lh --i ../surf/lh_DECIMATE_AREA_5.white --o ../surf/lh_DECIMATE_AREA_5.pial.T1 --pial --nsmooth 0 --rip-label ../label/lh_DECIMATE_AREA_5.cortex+hipamyg.label --pin-medial-wall ../label/lh_DECIMATE_AREA_5.cortex.label --aparc ../label/lh_DECIMATE_AREA_5.aparc.annot --repulse-surf ../surf/lh_DECIMATE_AREA_5.white --white-surf ../surf/lh_DECIMATE_AREA_5.white"

export SUBJECTS_DIR=$SUBJECTS_DIR/subject
test_command "cd $SUBJECTS_DIR/mri && echo "$cmd_1_banner" && $cmd_1 && echo "$cmd_2_banner" && $cmd_2 && echo "$cmd_3_banner" && $cmd_3"

export SUBJECTS_DIR=$SUBJECTS_DIR

# cmd 1 diff
(cd $SUBJECTS_DIR && ls -lt surf/lh_DECIMATE_AREA_5.white.preaparc* && ../../../mris_diff/mris_diff surf/lh_DECIMATE_AREA_5.white.preaparc surf/lh_DECIMATE_AREA_5.white.preaparc.REF_PLACE_SURFACES)

# cmd 2 diff
(cd $SUBJECTS_DIR && ls -lt surf/lh_DECIMATE_AREA_5.white* && ../../../mris_diff/mris_diff surf/lh_DECIMATE_AREA_5.white surf/lh_DECIMATE_AREA_5.white.REF_PLACE_SURFACES)

# cmd 3 diff
(cd $SUBJECTS_DIR && ls -lt surf/lh_DECIMATE_AREA_5.pial.T1* && ../../../mris_diff/mris_diff surf/lh_DECIMATE_AREA_5.pial.T1 surf/lh_DECIMATE_AREA_5.pial.T1.REF_PLACE_SURFACES)

