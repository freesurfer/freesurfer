#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$1"

cd tetgen

make && make tetlib

mkdir -p ${INSTALL_DIR}/lib ${INSTALL_DIR}/include ${INSTALL_DIR}/bin
cp -f tetgen ${INSTALL_DIR}/bin/
cp -f libtet.a ${INSTALL_DIR}/lib/
cp -f tetgen.h ${INSTALL_DIR}/include/
