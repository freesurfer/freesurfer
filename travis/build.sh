#!/usr/bin/env bash
set -ex

#
# script to run the travis build step
#

./setup_configure &>> build.log

config_flags="--with-pkgs-dir=${PWD}/build-packages --disable-Werror --disable-GUI-build"
[[ "$TRAVIS_OS_NAME" == "osx" ]] && config_flags="${configure_flags} F77=/usr/local/bin/gfortran-4.9 CC=/usr/local/bin/gcc-4.9 CXX=/usr/local/bin/g++-4.9"

./configure ${config_flags} &>> build.log
make -j4 &>> build.log
