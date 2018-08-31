#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

cd tetgen

make && make tetlib

mkdir -p ${INSTALL_DIR}/lib ${INSTALL_DIR}/include ${INSTALL_DIR}/bin
cp -f tetgen ${INSTALL_DIR}/bin/
cp -f libtet.a ${INSTALL_DIR}/lib/
cp -f tetgen.h ${INSTALL_DIR}/include/
