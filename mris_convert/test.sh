#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# don't remove test output before each test_command 
FSTEST_NO_DATA_RESET=1 && init_testdata

test_command mris_convert ref/lh.white ./lh.white.asc
test_command mris_convert ./lh.white.asc ./lh.new-white
test_command mris_convert ref/lh.white ./lh.white.gii
test_command mris_convert ./lh.white.gii ./lh.gifti-white
test_command mris_convert ref/lh.white ./lh.white.vtk
test_command mris_convert ./lh.white.vtk ./lh.vtk-white

test_command mris_convert -c ref/lh.thickness ref/lh.white ./lh.new-thickness.asc
test_command mris_convert -c ref/lh.thickness ref/lh.white ./lh.thickness.mgh
test_command mris_convert -c ./lh.thickness.mgh ref/lh.white ./lh.mgh-thickness.asc
test_command mris_convert -c ref/lh.thickness ref/lh.white ./lh.thickness.gii
test_command mris_convert -c ref/lh.thickness ref/lh.white ./lh.white.thickness.vtk
test_command mris_convert -c ref/lh.thickness.asc ref/lh.white ./lh.new-thickness
test_command mris_convert -c ./lh.new-thickness ref/lh.white ./lh.newer-thickness.asc
test_command mris_convert -c ./lh.thickness.gii ./lh.white.gii ./lh.gifti-thickness.asc
test_command mris_convert -c ./lh.white.thickness.vtk ./lh.white.vtk ./lh.vtk-thickness.asc

compare_surf ref/lh.white ./lh.new-white
compare_surf ref/lh.white ./lh.gifti-white
compare_surf ref/lh.white ./lh.vtk-white

compare_file ./lh.new-thickness.asc ./lh.newer-thickness.asc
compare_file ref/lh.thickness.asc ./lh.vtk-thickness.asc
compare_file ref/lh.thickness.asc ./lh.gifti-thickness.asc
compare_file ref/lh.thickness.asc ./lh.mgh-thickness.asc
compare_file ref/lh.thickness.asc ./lh.new-thickness.asc

test_command mris_convert ref/lh.caret.gifti.ascii.pial.gii ./lh.caret.gifti.ascii.pial
test_command mris_convert ref/lh.caret.gifti.base64.pial.gii ./lh.caret.gifti.base64.pial
test_command mris_convert ref/lh.caret.gifti.gzip_base64.pial.gii ./lh.caret.gifti.gzip_base64.pial
test_command mris_convert ref/lh.caret.gifti.fs.pial ./lh.caret.gifti.fs.pial.gii
test_command mris_convert lh.caret.gifti.fs.pial.gii ./lh.caret.gifti.fs.new-pial

test_command mris_convert -c ref/lh.caret.gifti.gzip_base64.shape.gii ref/lh.caret.gifti.base64.pial.gii ./lh.caret.gifti.fs.shape.asc
test_command mris_convert -c ref/rh.colin.func.gii ref/rh.colin.fudicial.gii ./rh.colin.func.mgh

compare_surf ref/lh.caret.gifti.fs.pial ./lh.caret.gifti.ascii.pial
compare_surf ./lh.caret.gifti.ascii.pial ./lh.caret.gifti.base64.pial
compare_surf ./lh.caret.gifti.base64.pial ./lh.caret.gifti.gzip_base64.pial
compare_surf ref/lh.caret.gifti.fs.pial ./lh.caret.gifti.fs.new-pial
compare_file ref/lh.caret.gifti.shape.expected.asc ./lh.caret.gifti.fs.shape.asc

#
#
# test case: convert .annot <=> GIFTI LabelTable+NIFTI_INTENT_LABEL
# test/ref data: 
#     bert/label/rh.mpm.vpnl.annot
#     bert/surf/rh.white
test_command mris_convert --annot bert/label/rh.mpm.vpnl.annot bert/surf/rh.white ./rh.mpm.vpnl.annot.gii
test_command mris_convert --annot ./rh.mpm.vpnl.annot.gii bert/surf/rh.white ./rh.mpm.vpnl.annot.gii.annot

# annotation and ctab should be the same
compare_annot ./rh.mpm.vpnl.annot.gii.annot bert/label/rh.mpm.vpnl.annot --diff-ctab

#
#
# test case: convert .annot <=> GIFTI NIFTI_INTENT_RGBA_VECTOR
test_command mris_convert --no-writect --annot bert/label/rh.mpm.vpnl.annot bert/surf/rh.white ./rh-rgba.mpm.vpnl.annot.gii
test_command mris_convert --annot ./rh-rgba.mpm.vpnl.annot.gii bert/surf/rh.white ./rh-rgba.mpm.vpnl.annot.gii.annot

# annotation should be the same, ctab differences are expected
# 0 Unknown_Label_0 0 (0 0 0 0)
# 1 Unknown_Label_1 255 (255 0 0 0)              hOc2   255     0     0     0
# 2 Unknown_Label_2 25600 (0 100 0 0)            hOc1     0   100     0     0
# 3 Unknown_Label_3 65535 (255 255 0 0)          hOc4v  255   255     0     0
# 4 Unknown_Label_4 1376057 (57 255 20 0)        FG1     57   255    20     0
# 5 Unknown_Label_5 1705837 (109 7 26 0)         FG4    109     7    26     0
# 6 Unknown_Label_6 16711680 (0 0 255 0)         FG3      0     0   255     0
# 7 Unknown_Label_7 16711935 (255 0 255 0)       FG2    255     0   255     0
# 8 Unknown_Label_8 16776960 (0 255 255 0)       hOc3v    0   255   255     0
compare_annot ./rh-rgba.mpm.vpnl.annot.gii.annot bert/label/rh.mpm.vpnl.annot

#
#
# test case: surface with group_avg_surface_area data <=> .gii
# test/ref data:
#     fsaverage/surf/rh.white
#     fsaverage/surf/rh-noextra.white
test_command mris_convert fsaverage/surf/rh.white ./rh-fsaverage.white.gii
test_command TRIANGULARSURFACE_NOEXTRA_WRITE=1 mris_convert ./rh-fsaverage.white.gii ./rh-fsaverage.noextra.white.gii.white

# compare rh-fsaverage.noextra* binaries w/o 'created by' and cmd lines
compare_bin ./rh-fsaverage.noextra.white.gii.white fsaverage/surf/rh-fsaverage.noextra.white

#
#
# test case: mris_convert --mergegifti/--splitgifti
# test/ref data:
#     bert/surf/rh-bert.noextra.white
#     bert/surf/rh.curv, bert/surf/rh.sulc, bert/surf/rh.thickness, bert/surf/rh.area
test_command mris_convert \
             --mergegifti \
             -c bert/surf/rh.curv bert/surf/rh.sulc bert/surf/rh.thickness bert/surf/rh.area \
             bert/surf/rh.white ./rh-bert.white+curv+sulc+thickness+area.gii
test_command mris_convert \
             --splitgifti --giftioutdir giftisplitout \
             ./rh-bert.white+curv+sulc+thickness+area.gii ./rh-bert.white.gii-merged.white
test_command TRIANGULARSURFACE_NOEXTRA_WRITE=1 mris_convert ./rh-bert.white.gii-merged.white ./rh-bert.noexta.white.gii-merged.white
test_command mris_convert -c giftisplitout/rh.curv.gii bert/surf/rh.white giftisplitout/rh.curv.gii.curv
test_command mris_convert -c giftisplitout/rh.sulc.gii bert/surf/rh.white giftisplitout/rh.sulc.gii.sulc
test_command mris_convert -c giftisplitout/rh.thickness.gii bert/surf/rh.white giftisplitout/rh.thickness.gii.thickness
test_command mris_convert -c giftisplitout/rh.area.gii bert/surf/rh.white giftisplitout/rh.area.gii.area

# compare rh-bert.noextra* binaries w/o 'created by' and cmd lines
compare_bin ./rh-bert.noexta.white.gii-merged.white bert/surf/rh-bert.noextra.white
# compare giftisplitout/rh.*.gii.[curv|sulc|thickness|area] with original bert/surf/rh.[curv|sulc|thickness|area]
compare_bin giftisplitout/rh.curv.gii.curv bert/surf/rh.curv
compare_bin giftisplitout/rh.sulc.gii.sulc bert/surf/rh.sulc
compare_bin giftisplitout/rh.thickness.gii.thickness bert/surf/rh.thickness
compare_bin giftisplitout/rh.area.gii.area bert/surf/rh.area
