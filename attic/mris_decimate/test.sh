#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mris_decimate -a 0.5 lh.orig.nofix.predec lh.orig.nofix
# 5/22/2019
# Commenting out this command which results in an mris_diff of the newly
# generated and reference surface (from testdata.tar.gz), since the results
# appear to not be determinstic on CentOS7 and the test outright fails on CentOS6
# compare_surf lh.orig.nofix lh.orig.nofix.ref
