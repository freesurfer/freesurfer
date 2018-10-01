#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

cd xml2
./configure --prefix=${INSTALL_DIR} --with-sax1 --disable-shared --with-minimum

make -j8
make install

