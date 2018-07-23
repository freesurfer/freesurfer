#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

cd netcdf/src

# CXX=null prevents building the cpp libraries (which we don't need)
./configure --prefix=${INSTALL_DIR} CXX=null

make  # can't be run in parallel
make install
