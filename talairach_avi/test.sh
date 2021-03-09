#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# this is a bit of a needy test, since it requires specific data to be in FREESURFER_HOME
# and it uses binaries and scripts scattered across the source (and build) tree

FSTEST_NO_DATA_RESET=1 && init_testdata

# talairach avi calls a few different binaries
for dir in io mri_convert mri_info mri_matrix_multiply mri_robust_register lta_convert; do
    export PATH="$(find_path $FSTEST_CWD $dir):$PATH"
done

# also calls items in the scripts directory
export PATH="$(find_path $FSTEST_SCRIPT_DIR scripts):$PATH"

# set up average link in the testing directory to mimic the FREESURFER_HOME
export FREESURFER_HOME="$FSTEST_TESTDATA_DIR"
ln -s $FSTEST_SCRIPT_DIR average

test_command talairach_avi --i nu.mgz --xfm tal.xfm --atlas 3T18yoSchwartzReactN32_as_orig

# convert .xfm to .lta since lta_diff no longer handles .xfm files
$(find_path $FSTEST_CWD lta_convert/lta_convert) --inmni tal.xfm --outlta tal.lta --src nu.mgz --trg mni305.cor.mgz

compare_lta tal.lta reference-tal.lta
