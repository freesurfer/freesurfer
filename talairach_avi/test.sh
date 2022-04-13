#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# this is a bit of a needy test, since it requires specific data to be in FREESURFER_HOME
# and it uses binaries and scripts scattered across the source (and build) tree

# don't reset the testdata directory after every command
FSTEST_NO_DATA_RESET=1 && init_testdata

# set up average link in the testing directory
ln -s $FSTEST_SCRIPT_DIR average

test_command talairach_avi --i nu.mgz --xfm tal.xfm --atlas 3T18yoSchwartzReactN32_as_orig

# convert .xfm to .lta since lta_diff no longer handles .xfm files
$(find_path $FSTEST_CWD lta_convert/lta_convert) --inmni tal.xfm --outlta tal.lta --src nu.mgz --trg mni305.cor.mgz
compare_lta tal.lta reference-tal.lta
