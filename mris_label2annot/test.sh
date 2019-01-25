#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_label2annot --s bert --h lh --ctab  bert/label/aparc.annot.ctab --a myaparc --l bert/label/lh.cortex.label --nhits nhits.mgh
compare_file bert/label/lh.myaparc.annot bert/label/lh.myaparc.annot.ref
