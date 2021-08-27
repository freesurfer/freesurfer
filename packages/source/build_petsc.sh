#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=${CC:-$(which gcc)}
export CXX=${CXX:-$(which g++)}

export CFLAGS=${CFLAGS:-"-msse2 -mfpmath=sse"}
export CXXFLAGS=${CXXFLAGS:-"-msse2 -mfpmath=sse"}

cd petsc

export PETSC_DIR=$(pwd)

python2 config/configure.py \
  --with-debugging=no \
  --download-f-blas-lapack=0 \
  --download-mpich=1 \
  --with-mpi=1 \
  --with-x=0 \
  --with-gnu-copyright-code=0 \
  --with-shared=0 \
  COPTFLAGS='-O3' \
  CXXOPTFLAGS='-O3' \
  FOPTFLAGS='-O3'

if [ "$(uname -s)" == "Linux" ]; then
  PETSC_ARCH="linux-gnu-c-opt"
elif [ "$(uname -s)" == "Darwin" ]; then
  PETSC_ARCH="darwin$(uname -r)-c-opt"
fi
export PETSC_ARCH

make -j 8

# build mpich
cd externalpackages/mpich2-1.0.5p4
./configure
make

# install everything
mkdir -p ${INSTALL_DIR}/bin ${INSTALL_DIR}/lib ${INSTALL_DIR}/include

cp -f ${PETSC_ARCH}/bin/* ${INSTALL_DIR}/
cp -f ${PETSC_ARCH}/lib/* ${INSTALL_DIR}/lib
cp -f ${PETSC_ARCH}/include/* ${INSTALL_DIR}/include
make clean

cd ../..
cp -f lib/${PETSC_ARCH}/* ${INSTALL_DIR}/lib
cp -f bmake/${PETSC_ARCH}/* ${INSTALL_DIR}/include
cp -rf include/* ${INSTALL_DIR}/include
make clean

ranlib ${INSTALL_DIR}/lib/lib*.a
