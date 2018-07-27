#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$1"

cd xml2
./configure --prefix=${INSTALL_DIR} --with-sax1 --disable-shared --with-minimum

make -j8
make install

