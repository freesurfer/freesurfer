#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# always run a multithreaded test, but generate single-threaded reference data
if [ "$FSTEST_REGENERATE" != true ]; then
    export OMP_NUM_THREADS=8
fi

if [ "$host_os" == "macos12" ]; then
   TESTDATA_SUFFIX=".clang13"
fi

test_command mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 subj1 lh
if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "ubuntu22" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "macos10" ]] || [[ "$host_os" == "macos12" ]]; then
   compare_surf subj1/surf/lh.orig subj1/surf/lh.orig.ref${TESTDATA_SUFFIX}
else
   compare_surf subj1/surf/lh.orig subj1/surf/lh.orig.ref
fi

test_command mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig -ga -seed 1234 subj3 rh
if [[ "$TESTDATA_SUFFIX" != "" ]] && [[ "$host_os" == "ubuntu20" ]] || [[ "$host_os" == "ubuntu22" ]] || [[ "$host_os" == "centos8" ]] || [[ "$host_os" == "macos10" ]] || [[ "$host_os" == "macos12" ]]; then
   compare_surf subj3/surf/rh.orig subj3/surf/rh.orig.ref${TESTDATA_SUFFIX}
else
   compare_surf subj3/surf/rh.orig subj3/surf/rh.orig.ref
fi

