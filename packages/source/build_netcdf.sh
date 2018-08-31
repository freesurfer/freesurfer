#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

cd netcdf/src

# CXX=null prevents building the cpp libraries (which we don't need)
./configure --prefix=${INSTALL_DIR} CXX=null

make  # can't be run in parallel
make install
