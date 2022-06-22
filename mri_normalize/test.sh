#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# mprage
test_command mri_normalize -mprage nu.mgz T1.mgz
compare_vol T1.mgz T1.ref.mgz

# aseg
test_command mri_normalize -aseg aseg.presurf.mgz norm.mgz brain.mgz

if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "macos10" ]]; then
   compare_vol brain.mgz brain.ref${TESTDATA_SUFFIX}.mgz
else
   compare_vol brain.mgz brain.ref.mgz
fi

# gentle
test_command mri_normalize -gentle nu.mgz gentle.mgz
compare_vol gentle.mgz gentle.ref.mgz

