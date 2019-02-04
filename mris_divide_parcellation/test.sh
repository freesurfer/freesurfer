#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command "mris_divide_parcellation bert lh aparc.annot splittable.txt split-sfg+pc 2>&1 | tee -a divide_parcel.out.raw && cat divide_parcel.out.raw | grep -v 'supposed to be reproducible' > divide_parcel.out"
compare_file divide_parcel.out divide_parcel.out.ref
