#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command samseg --i input.mgz --o output --threads 1 --options config.json
compare_vol output/seg.mgz seg.ref.mgz
