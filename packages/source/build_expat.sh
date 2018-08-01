#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

cd expat
./buildconf.sh && ./configure --prefix=${INSTALL_DIR}

make -j8
make install

