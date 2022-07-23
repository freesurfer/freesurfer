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

compare_vol bert/surf/lh.thickness bert/surf/lh.ref_thickness

