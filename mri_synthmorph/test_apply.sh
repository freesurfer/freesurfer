#!/usr/bin/env bash
. "$(dirname "$0")/../test.sh"

t() { test_command mri_synthmorph "$@" ; }

# image interpolation
t apply affine.lta moving.mgz out.mgz
compare_vol out.mgz affine.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3

# NIfTI format
t apply affine.lta moving.mgz out.nii.gz
compare_vol out.nii.gz affine.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3

# dense transform
t apply identity.mgz moving.mgz out.mgz -t uint8
compare_vol out.mgz moving.mgz --res-thresh 1e-3 --geo-thresh 1e-3

# matrix update
t apply -H rigid.lta -t uint8 moving.mgz out.mgz
compare_vol out.mgz header.mgz --res-thresh 1e-3 --geo-thresh 1e-3

# label interpolation
t apply -m nearest -t uint8 affine.lta labels.moving.mgz out.mgz
compare_vol out.mgz labels.affine.mgz --res-thresh 1e-3 --geo-thresh 1e-3

# multiple input pairs
t apply rigid.lta moving.mgz out_1.mgz moving.mgz out_2.mgz
compare_vol out_1.mgz rigid.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3
compare_vol out_2.mgz rigid.mgz --thresh 0.02 --res-thresh 1e-3 --geo-thresh 1e-3

# data types, command abbreviation
t a -t uint8 rigid.lta moving.mgz out.mgz
run_comparison mri_info --type out.mgz | grep -Fx uchar

t ap -t uint16 rigid.lta moving.mgz out.mgz
run_comparison mri_info --type out.mgz | grep -Fx ushrt

t app -t int16 rigid.lta moving.mgz out.mgz
run_comparison mri_info --type out.mgz | grep -Fx short

t appl -t int32 rigid.lta moving.mgz out.mgz
run_comparison mri_info --type out.mgz | grep -Fx int

t apply -t float32 rigid.lta moving.mgz out.mgz
run_comparison mri_info --type out.mgz | grep -Fx float

# method abbreviation
FSTEST_NO_DATA_RESET=1
for method in l li lin line linea linear n ne nea near neare neares nearest; do
    t apply -m "$method" rigid.lta moving.mgz out.mgz
done

# usage, help
t
t -h
t apply
t apply -h

# NIfTI warp
FSTEST_NO_DATA_RESET=1
mri_convert identity.mgz identity.nii.gz
t apply identity.nii.gz moving.mgz out.mgz -t uint8
compare_vol out.mgz moving.mgz --res-thresh 1e-3 --geo-thresh 1e-3

# illegal arguments
EXPECT_FAILURE=1
t slice
t apply -H identity.mgz moving.mgz fixed.mgz
t apply affine.lta moving.mgz out.mgz odd-number-of-io-pairs.mgz
