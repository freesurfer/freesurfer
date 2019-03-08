#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_smooth -n 1 -nw lh.orig.nofix lh.smoothwm.nofix
compare_surf lh.smoothwm.nofix lh.smoothwm.nofix.ref
