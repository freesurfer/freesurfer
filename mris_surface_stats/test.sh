#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command "mris_surface_stats -mask bert/label/lh.cortex.label -nsmooth 60 -surf_name bert/surf/lh.white bert/surf/lh.thickness | tee out.log"
compare_file out.log ref.log
