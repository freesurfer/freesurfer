#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_concat std.rh.*.mgh --o rhout.mgh
compare_vol rhout.mgh rhout.ref.mgh
