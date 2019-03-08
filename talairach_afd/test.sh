#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command "talairach_afd -xfm good_tal.xfm -afd ${FREESURFER_HOME}/fsafd -V | tee output"
compare_file output good_tal.txt -IfsafdDir

EXPECT_FAILURE=1 test_command "talairach_afd -xfm bad_tal.xfm -afd ${FREESURFER_HOME}/fsafd -V | tee output"
compare_file output bad_tal.txt -IfsafdDir
