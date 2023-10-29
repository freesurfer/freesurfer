#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_curvature_stats -o raw -h 10 -G cvs_avg35 rh rh.curv

if [ "$host_os" == "macos12" ]; then
   TESTDATA_SUFFIX=".clang13"
fi

compare_file raw ref${TESTDATA_SUFFIX}
for stat in BE C FI H K1 K2 K S SI; do
    if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] ||  [[ "$host_os" == "ubuntu22" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "centos9" ]] || [[ "$host_os" == "macos10" ]] || [[ "$host_os" == "macos12" ]]; then
       compare_file raw.rh.smoothwm.${stat}.hist ref${TESTDATA_SUFFIX}.rh.smoothwm.${stat}.hist
    else
       compare_file raw.rh.smoothwm.${stat}.hist ref.rh.smoothwm.${stat}.hist
    fi
done

