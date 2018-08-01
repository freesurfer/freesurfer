#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

export CC=$(which gcc)
export CXX=$(which g++)

# building minc requires netcdf, if NETCDF_DIR is not set, hope that it
# has been built right above the INSTALL_DIR path
if [ -z "${NETCDF_DIR}" ] ; then
  packages="$(dirname $(dirname ${INSTALL_DIR}))"
  # search for the header
  for header in $(find ${packages} -path "*/include/netcdf.h") ; do
    # check to make sure this potential netcdf directory has libnetcdf.a
    toplevel=$(dirname $(dirname $header))
    [ -e "${toplevel}/lib/libnetcdf.a" ] && NETCDF_DIR=${toplevel}
  done
fi

if [ ! -e "${NETCDF_DIR}" ] ; then
  echo "error: could not find netcdf installation - set NETCDF_DIR to a valid install path"
  exit 1
fi

cd minc-1.5

./configure --prefix=${INSTALL_DIR} CPPFLAGS=-I${NETCDF_DIR}/include LDFLAGS=-L${NETCDF_DIR}/lib

make -j8
make install
