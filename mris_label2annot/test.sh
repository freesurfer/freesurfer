#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

FSTEST_NO_DATA_RESET=1 && init_testdata

#
#
# test case: labels specified with --l are all found in colortab, and labels in colortab are sequentially numbered
# test/ref data:
#     bert/label/colortable_vpnl.txt
#     bert/label/rh.FG1.mpm.vpnl.label
#     bert/label/rh.FG2.mpm.vpnl.label
#     bert/label/rh.FG3.mpm.vpnl.label
#     bert/label/rh.FG4.mpm.vpnl.label
#     bert/label/rh.hOc1.mpm.vpnl.label
#     bert/label/rh.hOc2.mpm.vpnl.label
#     bert/label/rh.hOc3v.mpm.vpnl.label
#     bert/label/rh.hOc4v.mpm.vpnl.label
#     bert/label/rh.mpm.vpnl.ref.annot
test_command mris_label2annot \
                 --s bert --ctab bert/label/colortable_vpnl.txt --hemi rh --annot-path annotoutdir/rh.mpm.vpnl.test.annot --maxstatwinner --noverbose \
                 --l bert/label/rh.FG1.mpm.vpnl.label   \
                 --l bert/label/rh.FG2.mpm.vpnl.label   \
                 --l bert/label/rh.FG3.mpm.vpnl.label   \
                 --l bert/label/rh.FG4.mpm.vpnl.label   \
                 --l bert/label/rh.hOc1.mpm.vpnl.label  \
                 --l bert/label/rh.hOc2.mpm.vpnl.label  \
                 --l bert/label/rh.hOc3v.mpm.vpnl.label \
                 --l bert/label/rh.hOc4v.mpm.vpnl.label
compare_annot annotoutdir/rh.mpm.vpnl.test.annot bert/label/rh.mpm.vpnl.ref.annot


#
#
# test case: labels specified with --l are all found in colortab, and labels in colortab are not sequentially numbered
# test/ref data:
#     bert/label/colortable_vpnl.skippedlabelid.txt
#     bert/label/rh.FG1.mpm.vpnl.label
#     bert/label/rh.FG2.mpm.vpnl.label
#     bert/label/rh.FG3.mpm.vpnl.label
#     bert/label/rh.FG4.mpm.vpnl.label
#     bert/label/rh.hOc1.mpm.vpnl.label
#     bert/label/rh.hOc2.mpm.vpnl.label
#     bert/label/rh.hOc3v.mpm.vpnl.label
#     bert/label/rh.hOc4v.mpm.vpnl.label
#     bert/label/rh.mpm.vpnl.ref.annot
test_command mris_label2annot \
                 --s bert --ctab bert/label/colortable_vpnl.skippedlabelid.txt --hemi rh --annot-path annotoutdir/rh.mpm.vpnl.test.skippedlabelid.annot --maxstatwinner --noverbose \
                 --l bert/label/rh.FG1.mpm.vpnl.label   \
                 --l bert/label/rh.FG2.mpm.vpnl.label   \
                 --l bert/label/rh.FG3.mpm.vpnl.label   \
                 --l bert/label/rh.FG4.mpm.vpnl.label   \
                 --l bert/label/rh.hOc1.mpm.vpnl.label  \
                 --l bert/label/rh.hOc2.mpm.vpnl.label  \
                 --l bert/label/rh.hOc3v.mpm.vpnl.label \
                 --l bert/label/rh.hOc4v.mpm.vpnl.label
compare_annot annotoutdir/rh.mpm.vpnl.test.skippedlabelid.annot bert/label/rh.mpm.vpnl.ref.annot


#
#
# test case: --ldir <>, all labels are found in colortab, and labels in colortab are sequentially numbered
# test/ref data:
#     bert/label/colortable_vpnl.txt
#     bert/label/rh.FG1.mpm.vpnl.label
#     bert/label/rh.FG2.mpm.vpnl.label
#     bert/label/rh.FG3.mpm.vpnl.label
#     bert/label/rh.FG4.mpm.vpnl.label
#     bert/label/rh.hOc1.mpm.vpnl.label
#     bert/label/rh.hOc2.mpm.vpnl.label
#     bert/label/rh.hOc3v.mpm.vpnl.label
#     bert/label/rh.hOc4v.mpm.vpnl.label
#     bert/label/rh.mpm.vpnl.ref.annot
test_command mris_label2annot --s bert --ctab bert/label/colortable_vpnl.txt --hemi rh --annot-path annotoutdir/rh.mpm.vpnl.test.ldir.annot --maxstatwinner --noverbose --ldir bert/label
compare_annot annotoutdir/rh.mpm.vpnl.test.ldir.annot bert/label/rh.mpm.vpnl.ref.annot


