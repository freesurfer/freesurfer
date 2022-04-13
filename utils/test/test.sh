#!/usr/bin/env bash
source "$(dirname $0)/../../test.sh"

export PATH="$FSTEST_CWD:$PATH"

test_command testcolortab ${FREESURFER_HOME}/FreeSurferColorLUT.txt
test_command test_c_nr_wrapper
test_command extest
test_command inftest
test_command test_TriangleFile_readWrite
test_command topology_test
test_command tiff_write_image
test_command sc_test
test_command sse_mathfun_test
