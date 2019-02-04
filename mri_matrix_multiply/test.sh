#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_matrix_multiply -v -im bert/mri/transforms/talairach.xfm -im $SUBJECTS_DIR/cvs_avg35/mri/transforms/talairach.xfm -v -om mmult.xfm
compare_file mmult.xfm mmult.xfm.ref