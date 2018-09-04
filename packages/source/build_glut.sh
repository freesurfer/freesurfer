#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

if [ "$(uname -s)" == "Darwin" ] ; then
  export  CFLAGS="-I/usr/X11R6/include"
  export LDFLAGS="-L/usr/X11R6/lib"
fi

cd glut
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} .

make -j8
make install

