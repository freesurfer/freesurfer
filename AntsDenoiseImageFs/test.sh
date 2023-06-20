#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command ../../mri_convert/mri_convert -dsold 6 6 6 -i T1.mgz -o T1_downsample.mgz && AntsDenoiseImageFs -i T1_downsample.mgz -o T1.out.mgz
compare_vol T1.ref.mgz T1.out.mgz

