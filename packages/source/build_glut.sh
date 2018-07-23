#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

cd glut
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} .

make -j8
make install

