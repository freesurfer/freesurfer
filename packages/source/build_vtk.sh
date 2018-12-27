#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

cd vtk-v5.10.1

# CMAKE_PREFIX_PATH is needed to find the built tcltktixblt 
#   Without it, it is either not found, or the one in /usr/lib etc. may be the wrong version
#
cmake . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS="-DGLX_GLXEXT_LEGACY" \
  -DCMAKE_CXX_FLAGS="-DGLX_GLXEXT_LEGACY" \
  -DCMAKE_PREFIX_PATH=${INSTALL_DIR}/../../tcltktixblt/8.4.6/ \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DVTK_WRAP_TCL=ON

make -j8
make install
