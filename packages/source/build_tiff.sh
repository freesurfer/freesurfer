#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

cd tiff
./configure --prefix=${INSTALL_DIR} --noninteractive

make -j8
make install

