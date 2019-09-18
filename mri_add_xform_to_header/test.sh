#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# use talairach.xfm transform
mri_add_xform_to_header -v -c talairach.xfm orig.mgz out.mgz
compare_vol orig.mgz orig.ref.mgz

