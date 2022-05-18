#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_curvature_stats -o raw -h 10 -G cvs_avg35 rh rh.curv

compare_file raw ref${TESTDATA_SUFFIX}
for stat in BE C FI H K1 K2 K S SI; do
    compare_file raw.rh.smoothwm.${stat}.hist ref${TESTDATA_SUFFIX}.rh.smoothwm.${stat}.hist
done
