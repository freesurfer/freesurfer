#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

test_command samseg --i input.mgz --o output --threads 1 --options config.json
compare_vol output/seg.mgz seg.ref${TESTDATA_SUFFIX}.mgz
