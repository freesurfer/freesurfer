#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command "mri_segstats --subject bert --etiv-only | grep atlas_icv > etiv.txt"
compare_file etiv.txt etiv.ref.txt
