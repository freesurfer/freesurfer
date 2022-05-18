#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# mprage
test_command mri_normalize -mprage nu.mgz T1.mgz
compare_vol T1.mgz T1.ref.mgz

# aseg
test_command mri_normalize -aseg aseg.presurf.mgz norm.mgz brain.mgz
compare_vol brain.mgz brain.ref${TESTDATA_SUFFIX}.mgz

# gentle
test_command mri_normalize -gentle nu.mgz gentle.mgz
compare_vol gentle.mgz gentle.ref.mgz
