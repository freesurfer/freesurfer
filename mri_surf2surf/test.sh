#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_surf2surf \
    --hemi lh \
    --srcsubject bert \
    --srcsurfval thickness \
    --src_type curv \
    --trgsubject fsaverage \
    --trg_type curv \
    --trgsurfval bert/surf/lh.thickness

# compare_vol bert/surf/lh.thickness bert/surf/lh.ref_thickness
# 
# With reported error,
# ... Volumes differ in geometry row=1 col=1 diff=0.000000 (1.51904e-19)...
# Even a smaller threshold allows this to succeed (more than 18 leading zeroes to r.h.s of decimal point).
compare_vol --thresh .000000000000000000151904 bert/surf/lh.thickness bert/surf/lh.ref_thickness

