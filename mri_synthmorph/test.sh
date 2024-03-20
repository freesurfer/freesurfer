#!/usr/bin/env bash
. "$(dirname "$0")/../test.sh"

t() { test_command mri_synthmorph "$@" ; }

# affine inverse consistency
t -m affine -t out.lta -T inv.lta moving.mgz fixed.mgz
lta_diff out.lta inv.lta --invert2 | awk 'END {print $0; exit !($0<1e-3)}'

# rigid inverse consistency
t -m rigid -t out.lta -T inv.lta moving.mgz fixed.mgz
lta_diff out.lta inv.lta --invert2 | awk 'END {print $0; exit !($0<1e-3)}'

# rigid registration
t -m rigid -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz rigid.mgz --thresh 0.02

# geometry update
t -m rigid -Ho out.mgz moving.mgz fixed.mgz
compare_vol out.mgz header.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3

# affine registration
t -m affine -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz affine.mgz --thresh 0.02

# affine symmetry
t -m affine -O out.mgz fixed.mgz moving.mgz
compare_vol out.mgz affine.mgz --thresh 0.02

# deformable registration with initialization
t -m deform -i affine.lta -o out_1.mgz -O out_2.mgz moving.mgz fixed.mgz
compare_vol out_1.mgz deform_1.mgz --thresh 0.02
compare_vol out_2.mgz deform_2.mgz --thresh 0.02

# joint registration
t -m joint -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz joint.mgz --thresh 0.02

# default model
t -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz joint.mgz --thresh 0.02

# deformable flags
t moving.mgz fixed.mgz -mdeform -j16 -e256 -n7 -r0.5
t moving.mgz fixed.mgz -m joint -j 16 -e 192 -n 5 -r 0.7

# illegal arguments
EXPECT_FAILURE=1
t moving.mgz fixed.mgz -m banana
t moving.mgz fixed.mgz -e 1
t moving.mgz fixed.mgz -n 4
t moving.mgz fixed.mgz -r 0
t moving.mgz fixed.mgz -r 1
t moving.mgz fixed.mgz -Hm deform
t moving.mgz fixed.mgz -Hm joint
