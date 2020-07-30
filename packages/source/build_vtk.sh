#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

cd VTK-7.1.1

if [ "$(uname)" == "Darwin" ]; then
  EXTRA_OPTIONS="-DBUILD_SHARED_LIBS=OFF"
fi

cmake . $EXTRA_OPTIONS \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DVTK_RENDERING_BACKEND:STRING=OpenGL \
  -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}

cmake --build . --target all -j8
cmake --build . --target install
