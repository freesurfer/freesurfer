#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_parse_sdcmdir --sortbyrun --d . --o dicomdir.sumfile
compare_file dicomdir.sumfile dicomdir.ref.sumfile
