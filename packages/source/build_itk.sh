#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build_itk.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

mkdir build
cd build

cmake ../ITK-5.0a01 -G "Unix Makefiles" \
-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
-DITK_BUILD_DEFAULT_MODULES=OFF \
-DITKGroup_Core=ON \
-DITKGroup_Filtering=ON \
-DITKGroup_Segmentation=ON \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS="-msse2 -mfpmath=sse" \
-DCMAKE_C_FLAGS="-msse2 -mfpmath=sse"

make -j8
make install
