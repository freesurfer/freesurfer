#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_robust_template \
    --mov 001.mgz 001.mgz \
    --average 1 \
    --template rawavg.mgz \
    --satit \
    --inittp 1 \
    --fixtp \
    --noit \
    --iscale \
    --subsample 200 \

compare_vol rawavg.mgz 001.mgz --notallow-acq
