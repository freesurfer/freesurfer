#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_relabel_hypointensities mri/aseg.presurf.mgz surf aseg_hypoint_new.mgz
compare_vol aseg_hypoint_new.mgz aseg_hypoint_ref.mgz
