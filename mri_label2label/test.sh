#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_label2label \
    --srclabel bert/label/lh.BA1_exvivo.label \
    --srcsubject bert \
    --trglabel cvs_avg35/label/lh.BA1.label \
    --trgsubject cvs_avg35 \
    --regmethod surface \
    --hemi lh

compare_file cvs_avg35/label/lh.BA1.label label_ref/lh.BA1.label
