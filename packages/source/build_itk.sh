#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

mkdir build
cd build

cmake ../ITK \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-msse2 -mfpmath=sse" \
  -DCMAKE_C_FLAGS="-msse2 -mfpmath=sse" \
  -DITK_BUILD_DEFAULT_MODULES=OFF \
  -DITKGroup_Core=ON \
  -DITKGroup_Filtering=ON \
  -DITKGroup_Segmentation=ON \
  -DITKGroup_IO=ON \
  -DModule_AnisotropicDiffusionLBR=ON \
  -DBUILD_TESTING=OFF

cmake --build . --target all -j8
cmake --build . --target install
