#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_segment mri/brainmask.mgz wm_segment_out.mgz
compare_vol wm_segment_out.mgz wm_segment_ref.mgz
