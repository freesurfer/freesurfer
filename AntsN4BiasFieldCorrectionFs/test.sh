#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command ../../mri_convert/mri_convert -ds 6 6 6 -i T1.mgz -o T1_downsample.mgz && AntsN4BiasFieldCorrectionFs -i T1_downsample.mgz -o T1.out.mgz
if [ "$host_os" == "ubuntu18" ]; then
   compare_vol --thresh 0.00042725 T1.ref${TESTDATA_SUFFIX}.mgz T1.out.mgz
else
   compare_vol T1.ref${TESTDATA_SUFFIX}.mgz T1.out.mgz
fi
