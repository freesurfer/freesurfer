#!/usr/bin/env bash

# set -x

# switch to hardcoded values to run this script manually
# in sandbox without typing "make"

# target_dir=${PWD}/mri_synthmorph
target_dir=$1

# arch_prefix="mu40-crazy"
arch_prefix=$2

# (cd $target_dir && git annex unlock ${arch_prefix}.tgz)
(cd $target_dir && git annex unlock ${arch_prefix}.tgz > /dev/null 2>&1)
(cd $target_dir && rm -rf $arch_prefix && tar zxpf ${arch_prefix}.tgz)
(cd $target_dir && ls $arch_prefix > /dev/null 2>&1)
