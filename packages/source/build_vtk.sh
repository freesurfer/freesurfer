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

# NOTE: ITK 4.13 and uses cmake featues from cmake 2.12
# this older cmake < 3.3 version does not support
# respecting the CMAKE_CXX_STANDARD values,
# and therefore the compiler is self identified
# as C++14 for gcc 7.5.0+ compilers
#
# If using a compiler that defaults to C++14 or C++17
# then freesurfer proper will need to be compiled
# with at least that level of support
# -DCMAKE_CXX_STANDARD:STRING=14 or 17 as dictated
# by the compiler.
cmake ../VTK-7.1.1 $EXTRA_OPTIONS \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
  -DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
  -DCMAKE_CXX_STANDARD:STRING=11 \
  -DCMAKE_C_STANDARD:STRING=11 \
  -DVTK_RENDERING_BACKEND:STRING=OpenGL

cmake --build . --target all -j 8
cmake --build . --target install
