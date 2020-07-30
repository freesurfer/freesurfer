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
  -DCMAKE_CXX_STANDARD:STRING=11 \
  -DCMAKE_C_STANDARD:STRING=11 \
  -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
  -DITKGroup_Core:BOOL=ON \
  -DITKGroup_Filtering:BOOL=ON \
  -DITKGroup_Segmentation:BOOL=ON \
  -DITKGroup_IO:BOOL=ON \
  -DModule_AnisotropicDiffusionLBR:BOOL=ON \
  -DBUILD_TESTING:BOOL=OFF

cmake --build . --target all -j8
cmake --build . --target install
