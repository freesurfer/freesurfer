#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# conform
test_command mri_convert rawavg.mgz orig.mgz --conform
compare_vol orig.mgz orig.ref.mgz

# dicom
test_command mri_convert dcm/261000-10-60.dcm dicom.mgz
compare_vol dicom.mgz freesurfer.mgz

# nifti
test_command mri_convert nifti.nii nifti.mgz
compare_vol nifti.mgz freesurfer.mgz --notallow-acq --geo-thresh 0.000008

# analyze
test_command mri_convert analyze.img analyze.mgz
compare_vol analyze.mgz freesurfer.mgz --notallow-acq --geo-thresh 0.000008

# mri_make_uchar
test_command mri_make_uchar nu.mgz talairach.xfm nuuc.mgz
compare_vol nuuc.mgz nuuc.ref.mgz

# apply downsampled morph with atlas geometry that has odd number of slices
# (i.e. gcam->depth * gcam->spacing = gcam->atlas.depth - 1)
test_command mri_convert -at odd.m3z orig.mgz morphed.mgz
compare_vol morphed.mgz odd.ref.mgz
