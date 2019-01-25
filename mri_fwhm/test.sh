#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_fwhm --i HelixTensors.nii.gz --nframesmin 9 --auto-mask .2 --dat fwhm.dat
compare_file fwhm.dat fwhm_ref.dat
