#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi 
INSTALL_DIR="$(realpath $1)"

cd jpeg
./configure --prefix=${INSTALL_DIR}

make -j8
mkdir -p ${INSTALL_DIR}/lib ${INSTALL_DIR}/include && make install-lib

