#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# since we're comparing annotations, let's do the diff manually
# until there's an easier way
mris_diff=$(find_path $FSTEST_CWD mris_diff/mris_diff)

# create desikan parcellation
test_command mris_ca_label \
    -l bert/label/lh.cortex.label \
    -aseg bert/mri/aseg.mgz \
    -seed 1234 \
    bert lh bert/surf/lh.sphere.reg \
    ${FREESURFER_HOME}/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs \
    bert/label/lh.aparc.annot

$mris_diff --maxerrs 1000 --s1 bert --s2 bert --hemi lh \
    --aparc aparc --aparc2 aparc.reference

# create destrieux parcellation
test_command mris_ca_label \
    -l bert/label/lh.cortex.label \
    -aseg bert/mri/aseg.mgz \
    -seed 1234 \
    bert lh bert/surf/lh.sphere.reg \
    ${FREESURFER_HOME}/average/lh.destrieux.simple.2009-07-29.gcs \
    bert/label/lh.aparc.a2009s.annot

$mris_diff --maxerrs 1000 --s1 bert --s2 bert --hemi lh \
    --aparc aparc.a2009s --aparc2 aparc.a2009s.reference
