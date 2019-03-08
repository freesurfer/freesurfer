#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_annotation2label --subject bert --sd . --hemi rh --outdir labels

for label in ${FSTEST_TESTDATA_DIR}/labels/*.label; do
    compare_file $label labels_ref/$(basename $label)
done
