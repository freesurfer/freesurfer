#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_diff lh.sphere lh.sphere
EXPECT_FAILURE=1 test_command mris_diff lh.sphere rh.sphere
