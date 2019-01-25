#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_expand -v 111706 -thickness lh.white 0.5 lh.midgray
compare_surf lh.midgray lh.expected-midgray
