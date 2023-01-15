#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

if [ "$host_os" == "macos12" ]; then
   test_command mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz talairach.m3z \
    ${FSTEST_SCRIPT_DIR}/testdata/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.out.mgz
else
   test_command mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz talairach.m3z \
   ${FREESURFER_HOME}/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.out.mgz
fi

# The tests for mri_ca_label and AntsN4BiasFieldCorrectionFs were the only failures on
# MacOS 10.12 (Monterey) using reference testdata files with the .clang12 TESTDATA_SUFFIX
# - created for running tests on  MacOS10.15 (Catalina).  Given this and that MacOS 10.12
# uses clang13 (instead of clang12), then just directly amend the ifdef for each of these
# tests to look for newly generated reference files on MacOS 10.12 with .clang13 suffix
# So as of this writing TESTDATA_SUFFIX not defined for MacOS 12 and hardcode .clang13 below.
if [ "$host_os" == "macos12" ]; then
   compare_vol aseg.auto_noCCseg.out.mgz aseg.auto_noCCseg.clang13.mgz
elif [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "macos10" ]]; then
   compare_vol aseg.auto_noCCseg.out.mgz aseg.auto_noCCseg${TESTDATA_SUFFIX}.mgz
else
   compare_vol aseg.auto_noCCseg.out.mgz aseg.auto_noCCseg.mgz
fi

