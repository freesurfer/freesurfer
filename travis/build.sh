#!/usr/bin/env bash

#
# script to run the travis build step
# this lets us handle build failures better
#

function runstep {
  $@ >> build.log 2>&1 || { tail -n 50 build.log; echo "travis build failed during $@"; exit 1; }
}

runstep "./setup_configure"

config_flags="--with-pkgs-dir=${PWD}/build-packages --disable-Werror --disable-GUI-build"
[[ "$TRAVIS_OS_NAME" == "osx" ]] && config_flags="${config_flags} F77=/usr/local/bin/gfortran-4.9 CC=/usr/local/bin/gcc-4.9 CXX=/usr/local/bin/g++-4.9"

runstep ./configure ${config_flags}
runstep make -j4
