#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# don't remove test output before each test_command 
FSTEST_NO_DATA_RESET=1 && init_testdata

test_command mris_convert lh.white lh.white.asc
test_command mris_convert lh.white.asc lh.new-white
test_command mris_convert lh.white lh.white.gii
test_command mris_convert lh.white.gii lh.gifti-white
test_command mris_convert lh.white lh.white.vtk
test_command mris_convert lh.white.vtk lh.vtk-white

test_command mris_convert -c lh.thickness lh.white lh.new-thickness.asc
test_command mris_convert -c lh.thickness lh.white lh.thickness.mgh
test_command mris_convert -c lh.thickness.mgh lh.white lh.mgh-thickness.asc
test_command mris_convert -c lh.thickness lh.white lh.thickness.gii
test_command mris_convert -c lh.thickness lh.white lh.white.thickness.vtk
test_command mris_convert -c lh.thickness.asc lh.white lh.new-thickness
test_command mris_convert -c lh.new-thickness lh.white lh.newer-thickness.asc
test_command mris_convert -c lh.thickness.gii lh.white.gii lh.gifti-thickness.asc
test_command mris_convert -c lh.white.thickness.vtk lh.white.vtk lh.vtk-thickness.asc

compare_surf lh.white lh.new-white
compare_surf lh.white lh.gifti-white
compare_surf lh.white lh.vtk-white

compare_file lh.new-thickness.asc lh.newer-thickness.asc
compare_file lh.thickness.asc lh.vtk-thickness.asc
compare_file lh.thickness.asc lh.gifti-thickness.asc
compare_file lh.thickness.asc lh.mgh-thickness.asc
compare_file lh.thickness.asc lh.new-thickness.asc

mris_convert lh.caret.gifti.ascii.pial.gii lh.caret.gifti.ascii.pial
mris_convert lh.caret.gifti.base64.pial.gii lh.caret.gifti.base64.pial
mris_convert lh.caret.gifti.gzip_base64.pial.gii lh.caret.gifti.gzip_base64.pial
mris_convert lh.caret.gifti.fs.pial lh.caret.gifti.fs.pial.gii
mris_convert lh.caret.gifti.fs.pial.gii lh.caret.gifti.fs.new-pial

mris_convert -c lh.caret.gifti.gzip_base64.shape.gii lh.caret.gifti.base64.pial.gii lh.caret.gifti.fs.shape.asc
mris_convert -c rh.colin.func.gii rh.colin.fudicial.gii rh.colin.func.mgh

compare_surf lh.caret.gifti.fs.pial lh.caret.gifti.ascii.pial
compare_surf lh.caret.gifti.ascii.pial lh.caret.gifti.base64.pial
compare_surf lh.caret.gifti.base64.pial lh.caret.gifti.gzip_base64.pial
compare_surf lh.caret.gifti.fs.pial lh.caret.gifti.fs.new-pial
compare_file lh.caret.gifti.shape.expected.asc lh.caret.gifti.fs.shape.asc
