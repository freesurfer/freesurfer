#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_inflate rh.smoothwm.nofix rh.inflated.nofix
compare_surf rh.inflated.nofix rh.inflated.nofix.ref
