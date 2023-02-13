#!/usr/bin/env bash -x

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=${CC:-$(which gcc)}
export CXX=${CXX:-$(which g++)}

export CFLAGS=${CFLAGS:-"-msse2 -mfpmath=sse"}
export CXXFLAGS=${CXXFLAGS:-"-msse2 -mfpmath=sse"}

mkdir build
cd build

# NOTE: ITK 4.13 and uses cmake featues from cmake 2.12
# this older cmake < 3.3 version does not support
# respecting the CMAKE_CXX_STANDARD values,
# and therefore the compiler is self identified
# as C++14 for gcc 7.5.0+ compilers.
#
# If using a compiler that defaults to C++14 or C++17
# then freesurfer proper will need to be compiled
# with at least that level of support
# -DCMAKE_CXX_STANDARD:STRING=14 or 17 as dictated
# by the compiler.
#
# 2/2023 - Upgrading from 4.13.0 to 4.13.2, adding modules
# ITKGroup_Registration
# ITKGroup_Numerics 
# As of this writing, the top level ../../CMakeLists.txt still defines
# set(CMAKE_CXX_STANDARD 11)
# (see note above)
cmake ../ITK \
  -DCMAKE_VERBOSE_MAKEFILE=ON \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_CXX_FLAGS:STRING="${CXXFLAGS}" \
  -DCMAKE_C_FLAGS:STRING="${CFLAGS}" \
  -DCMAKE_CXX_STANDARD:STRING=11 \
  -DCMAKE_C_STANDARD:STRING=11 \
  -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF \
  -DITKGroup_Core:BOOL=ON \
  -DITKGroup_Filtering:BOOL=ON \
  -DITKGroup_Numerics:BOOL=ON \
  -DITKGroup_Registration:BOOL=ON \
  -DITKGroup_Segmentation:BOOL=ON \
  -DITKGroup_IO:BOOL=ON \
  -DModule_AnisotropicDiffusionLBR:BOOL=ON \
  -DBUILD_TESTING:BOOL=OFF

# cmake --build . --target all -j 8
# Preserve log output order, maybe easier to check for errors
cmake --build . --target all -j1
cmake --build . --target install
