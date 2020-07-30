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

if [ "$(uname)" == "Darwin" ]; then
  EXTRA_OPTIONS="-DBUILD_SHARED_LIBS=OFF"
fi

cmake ../VTK-7.1.1 $EXTRA_OPTIONS \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
  -DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
  -DCMAKE_CXX_STANDARD:STRING=11 \
  -DCMAKE_C_STANDARD:STRING=11 \
  -DVTK_RENDERING_BACKEND:STRING=OpenGL

cmake --build . --target all -j8
cmake --build . --target install
