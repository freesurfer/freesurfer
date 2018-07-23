#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build_vtk.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

cd vtk-v5.10.1

cmake . -DCMAKE_BUILD_TYPE=Release \
        -DVTK_WRAP_TCL=ON \
        -DCMAKE_C_FLAGS=-DGLX_GLXEXT_LEGACY \
        -DCMAKE_CXX_FLAGS=-DGLX_GLXEXT_LEGACY \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}

make -j8
make install
