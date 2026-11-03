#!/usr/bin/env bash
source "$(dirname "$0")/../test.sh"

# don't remove test output before each test_command 
FSTEST_NO_DATA_RESET=1 && init_testdata

t() { FSTEST_NO_DATA_RESET=1 && test_command move_pointset "$@" ; }

# $FSTEST_TESTDATA_DIR/move_pointset directory contents:
# TEST 1:
#    fix-to-mov.nii.gz      - Test 1 input, Synthmorph backward deformation fields
#    grid_mov.label         - Test 1 input, mov CRS grid
#    ref/ref_mov2fix.label  - Test 1 output reference, mov -> fix vox2vox mapping, generated from fix-to-mov.nii.gz

### TEST 1: test with Synthmorph backward deformation fields
# sample the deformation field @mov grid to warp coords from mov -> fix
t --labelfile move_pointset/grid_mov.label --warp move_pointset/fix-to-mov.nii.gz --outcoordspace voxel --o move_pointset/mov2fix_warp.label
compare_file move_pointset/mov2fix_warp.label move_pointset/ref/ref_mov2fix.label


