#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

if [ "$host_os" == "macos12" ]; then
   test_command mri_watershed -T1 -brain_atlas ${FSTEST_SCRIPT_DIR}/testdata/average/RB_all_withskull_2016-05-10.vc700.gca talairach_with_skull.lta T1.mgz brainmask.mgz
else
   test_command mri_watershed -T1 -brain_atlas ${FREESURFER_HOME}/average/RB_all_withskull_2016-05-10.vc700.gca talairach_with_skull.lta T1.mgz brainmask.mgz
fi

# Have not yet set TESSTDATA_SUFFIX as .clang13 for MacOS 12
if [ "$host_os" == "macos12" ]; then
   compare_vol brainmask.mgz brainmask.ref.clang13.mgz
elif [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "ubuntu22" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "centos9" ]] || [[ "$host_os" == "macos10" ]]; then
   compare_vol brainmask.mgz brainmask.ref${TESTDATA_SUFFIX}.mgz
else
   compare_vol brainmask.mgz brainmask.ref.mgz
fi

