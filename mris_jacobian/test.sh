#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_jacobian rh.white rh.sphere.reg rh.jacobian_white
compare_vol rh.jacobian_white rh.jacobian_white.ref
