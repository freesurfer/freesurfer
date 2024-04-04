#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

EXPECT_FAILURE=1 test_command "mri_diff lh.ribbon.mgz rh.ribbon.mgz --ssd --count | grep -v 'mri_diff' | tee diff.log"
compare_file diff.log diff_ref.log
