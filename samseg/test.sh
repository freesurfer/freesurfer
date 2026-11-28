#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

if [ "$host_os" == "macos12" ]; then
   if [ "$CLANG_VERSION" == "14" ]; then
      TESTDATA_SUFFIX=".clang14"
   else
      TESTDATA_SUFFIX=".clang13"
   fi
   # ensure TERM defined to known to be good value
   export TERM=linux
fi

test_command samseg --i input.mgz --o output --threads 1 --options config.json

if [[ "$host_os" == "ubuntu24" ]] || [[ "$host_os" == "centos9" ]] || [[ "$host_os" == "rocky9" ]]; then
      TESTDATA_SUFFIX=".gcc11"
fi

if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "ubuntu22" ]] || [[ "$host_os" == "ubuntu24" ]]; then
   compare_vol output/seg.mgz seg.ref${TESTDATA_SUFFIX}.ubuntu.mgz
elif [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "centos9" ]] || [[ "$host_os" == "rocky9" ]]; then
   compare_vol output/seg.mgz seg.ref${TESTDATA_SUFFIX}.centos.mgz
elif [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "macos10" ]] || [[ "$host_os" == "macos12" ]] ; then
   compare_vol output/seg.mgz seg.ref${TESTDATA_SUFFIX}.mgz
else
   compare_vol output/seg.mgz seg.ref.mgz
fi
