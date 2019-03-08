#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# relabel corpus callosum
test_command mri_fuse_segmentations \
    --aseg aseg.1.mgz \
    --nocc aseg.auto_noCCseg.1.mgz \
    --norm norm.1.mgz \
    norm.1.mgz cc.mgz
compare_vol cc.mgz cc.ref.mgz

# fuse segmentations using LTAs
test_command mri_fuse_segmentations \
    -a aseg.{1,2,3}.mgz \
    -c aseg.auto_noCCseg.{1,2,3}.mgz \
    -n norm.{1,2,3}.mgz \
    -t im{1,2,3}_to_base.lta \
    norm.base.mgz lta.mgz
compare_vol lta.mgz lta.ref.mgz

# use different sigma
test_command mri_fuse_segmentations \
    -a aseg.{1,2,3}.mgz \
    -c aseg.auto_noCCseg.{1,2,3}.mgz \
    -n norm.{1,2,3}.mgz \
    -t im{1,2,3}_to_base.lta \
    -s 4.0 \
    norm.base.mgz sigma.mgz
compare_vol sigma.mgz sigma.ref.mgz

# fuse segmentations using GCAMs and identity transform
test_command mri_fuse_segmentations \
    -a aseg.{1,2,3}.mgz \
    -c aseg.auto_noCCseg.{1,2,3}.mgz \
    -n norm.{1,2,3}.mgz \
    -t im1_to_im2.m3z identity.nofile im3_to_im2.m3z \
    norm.2.mgz gcam.mgz
compare_vol gcam.mgz gcam.ref.mgz
