#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_surf2vol --o vol.mgz --subject cvs_avg35 --so lh.white lh.thickness --so rh.white rh.thickness
compare_vol vol.mgz vol_ref.mgz