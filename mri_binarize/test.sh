#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_binarize --i aseg.mgz --o mask.mgz --match 17
compare_vol mask.mgz mask.ref.mgz
