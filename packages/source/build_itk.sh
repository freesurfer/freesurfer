#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=${CC:-$(which gcc)}
export CXX=${CXX:-$(which g++)}

export CFLAGS=${CFLAGS:-"-msse2 -mfpmath=sse"}
export CXXFLAGS=${CXXFLAGS:-"-msse2 -mfpmath=sse"}

mkdir build
cd build





cmake ../ITK \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
  -DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
  -DITK_BUILD_DEFAULT_MODULES=OFF \
  -DITKGroup_Core=ON \
  -DITKGroup_Filtering=ON \
  -DITKGroup_Segmentation=ON \
  -DITKGroup_IO=ON \
  -DModule_AnisotropicDiffusionLBR=ON \
  -DBUILD_TESTING=OFF

cmake --build . --target all -j8
cmake --build . --target install
