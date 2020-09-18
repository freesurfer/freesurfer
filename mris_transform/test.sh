#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# untaring is time consuming, so don't remove output before each test
FSTEST_NO_DATA_RESET=1 && init_testdata

# transform surface with LTA
test_command mris_transform lh.1.pial 1_to_2.lta lh.out.pial
compare_surf lh.out.pial lh.ref.pial

# transform surface with inverse LTA
test_command mris_transform --is-inverse lh.1.pial 2_to_1.lta lh.out.pial
compare_surf lh.out.pial lh.ref.pial

# transform surface with FSLMAT
test_command mris_transform --trx-src nu.1.mgz --trx-dst nu.2.mgz lh.1.pial 1_to_2.fslmat lh.out.pial
compare_surf lh.out.pial lh.ref.pial

# transform surface with FSLMAT
test_command mris_transform --trx-src nu.2.mgz --trx-dst nu.1.mgz --is-inverse lh.1.pial 2_to_1.fslmat lh.out.pial
compare_surf lh.out.pial lh.ref.pial

# transform surface with GCAM
test_command mris_transform --trx-dst nu.2.mgz lh.1.pial 1_to_2.m3z lh.out.pial
compare_surf lh.out.pial lh.ref.gcam.pial

# transform surface with inverse GCAM
test_command mris_transform --trx-dst nu.1.mgz --is-inverse lh.1.pial 2_to_1.m3z lh.out.pial
compare_surf lh.out.pial lh.ref.gcam.inv.pial

# just copy surface
test_command mris_transform lh.1.pial identity.nofile lh.out.pial
compare_surf lh.out.pial lh.1.pial
