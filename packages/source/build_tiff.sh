#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

export DIR_LIBGZ=/usr/lib  # make sure to compile with ZIP_SUPPORT

cd tiff
./configure --prefix=${INSTALL_DIR} --noninteractive

make -j8
make install

