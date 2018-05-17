#!/usr/bin/env bash
set -ex

#
# script to run the travis build step
#

./setup_configure

config_flags="--with-pkgs-dir=${PWD}/build-packages --disable-Werror --disable-GUI-build"
[[ "$TRAVIS_OS_NAME" == "osx" ]] && config_flags="${configure_flags} F77=/usr/local/bin/gfortran-4.9 CC=/usr/local/bin/gcc-4.9 CXX=/usr/local/bin/g++-4.9"

./configure ${config_flags}
make -j4
