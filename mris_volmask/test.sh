#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_volmask \
    --label_left_white 2 \
    --label_left_ribbon 3 \
    --label_right_white 41 \
    --label_right_ribbon 42 \
    --save_ribbon \
    bert

compare_vol bert/mri/ribbon.mgz bert/mri/ribbon.ref.mgz
