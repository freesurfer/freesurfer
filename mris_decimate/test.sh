#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_decimate -a 0.5 lh.orig.nofix.predec lh.orig.nofix
compare_surf lh.orig.nofix lh.orig.nofix.ref
