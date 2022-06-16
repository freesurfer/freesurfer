#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

test_command samseg --i input.mgz --o output --threads 1 --options config.json
if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu18" ]] || [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "macos10" ]]; then
   compare_vol output/seg.mgz seg.ref${TESTDATA_SUFFIX}.mgz
else
   compare_vol output/seg.mgz seg.ref.mgz
fi

