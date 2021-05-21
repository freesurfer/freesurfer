#!/usr/bin/env bash

set -e

# gts was added becuase it is used by the mri_decimate command

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=${CC:-$(which gcc)}
export CXX=${CXX:-$(which g++)}

export CFLAGS=${CFLAGS:-"-msse2 -mfpmath=sse"}
export CXXFLAGS=${CXXFLAGS:-"-msse2 -mfpmath=sse"}

cd gts-0.7.6

# Mac requires netpbm from homebrew
if [ "$(uname)" == "Darwin" ]; then
   ./configure --prefix=${INSTALL_DIR} CFLAGS="-DUSE_SURFACE_BTREE -I/usr/local/Cellar/netpbm/10.73.24/include/netpbm"
elif [ "$(uname)" == "Linux" ]; then
   ./configure --prefix=${INSTALL_DIR} CFLAGS=-DUSE_SURFACE_BTREE
fi

# build
make

# install
mkdir -p ${INSTALL_DIR}
echo "INSTALL_DIR=$INSTALL_DIR"
make install

