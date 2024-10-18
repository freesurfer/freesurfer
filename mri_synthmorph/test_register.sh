#!/usr/bin/env bash
. "$(dirname "$0")/../test.sh"

t() { test_command mri_synthmorph "$@" ; }

# affine registration
t -m affine -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz affine.mgz --thresh 0.02

# affine symmetry
t -m aff -O out.mgz fixed.mgz moving.mgz
compare_vol out.mgz affine.mgz --thresh 0.02

# affine inverse consistency
t -m a -t out.lta -T inv.lta moving.mgz fixed.mgz
lta_diff out.lta inv.lta --invert2 | awk 'END {print $0; exit !($0<1e-3)}'

# rigid inverse consistency
t -m rigid -t out.lta -T inv.lta moving.mgz fixed.mgz
lta_diff out.lta inv.lta --invert2 | awk 'END {print $0; exit !($0<1e-3)}'

# rigid registration
t -m rig -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz rigid.mgz --thresh 0.02

# geometry update
t -m r -Ho out.mgz moving.mgz fixed.mgz
compare_vol out.mgz header.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3

# deformable registration with initialization
t -m deform -i affine.lta -o out_1.mgz -O out_2.mgz moving.mgz fixed.mgz
compare_vol out_1.mgz deform_1.mgz --thresh 0.02
compare_vol out_2.mgz deform_2.mgz --thresh 0.02

# deformable registration with mid-space initialization
t -m def -Mi affine.lta -o out.nii.gz moving.mgz fixed.mgz
compare_vol out.nii.gz deform_mid.nii.gz --thresh 0.02

# joint registration
t -m joint -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz joint.mgz --thresh 0.02

# default model
t -o out.mgz moving.mgz fixed.mgz
compare_vol out.mgz joint.mgz --thresh 0.02

# help, usage, command abbreviation
FSTEST_NO_DATA_RESET=1
t register -h
for cmd in r re reg regi regis registe register; do
    t "$cmd"
done

# deformable flags, explicit command
t register moving.mgz fixed.mgz -md -j16 -e256 -n7 -r0.5
t register moving.mgz fixed.mgz -m j -j 16 -e 192 -n 5 -r 0.7

# NIfTI warps
FSTEST_NO_DATA_RESET=1
mri_convert=$(find_path $FSTEST_CWD mri_convert/mri_convert)
t -t out.nii.gz moving.mgz fixed.mgz
test_command $mri_convert -odt float -at out.nii.gz moving.mgz out.mgz
compare_vol out.mgz joint.mgz --thresh 1

# illegal arguments
EXPECT_FAILURE=1
t moving.mgz fixed.mgz -m banana
t moving.mgz fixed.mgz -e 1
t moving.mgz fixed.mgz -n 4
t moving.mgz fixed.mgz -r 0
t moving.mgz fixed.mgz -r 1
t moving.mgz fixed.mgz -Hm deform
t moving.mgz fixed.mgz -Hm joint
