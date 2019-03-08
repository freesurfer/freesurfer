#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_concatenate_lta -invert1 nick13_to_nickbase.lta identity.nofile nickbase_to_nick13.lta
compare_file nickbase_to_nick13.lta reference-nickbase_to_nick13.lta

test_command mri_concatenate_lta 002.lta tp_to_base.lta 002-long.lta
compare_file 002-long.lta reference-002-long.lta
