#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_anatomical_stats -log lh.stats bert lh
compare_file lh.stats lh.ref.stats

test_command mris_anatomical_stats -a aparc.annot -f stats.table bert lh
compare_file stats.table stats.ref.table -I# -bw

test_command mris_anatomical_stats -cortex bert/label/lh.cortex.label -crosscheck -a bert/label/lh.aparc.annot -f aparc.stats.table bert lh
compare_file aparc.stats.table aparc.stats.ref.table -I# -bw

test_command mris_anatomical_stats -cortex bert/label/lh.cortex.label -crosscheck -a bert/label/lh.aparc.a2009s.annot -f a2009s.stats.table bert lh
compare_file a2009s.stats.table a2009s.stats.ref.table -I# -bw
