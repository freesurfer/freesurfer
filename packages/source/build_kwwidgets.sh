#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

cd KWWidgets-HEAD-cvs

# building kwwidgets requires vtk, if VTK_DIR is not set, hope that it
# has been built right above the INSTALL_DIR path
if [ -z "${VTK_DIR}" ] ; then
  packages="$(dirname $(dirname ${INSTALL_DIR}))"
  # search for the header
  for cmakepath in $(find ${packages} -name "VTKConfig.cmake") ; do
    VTK_DIR="$(dirname ${cmakepath})"
  done
fi

if [ ! -e "${VTK_DIR}" ] ; then
  echo "error: could not find vtk installation - set VTK_DIR to a valid install path"
  exit 1
fi

cmake_cxx_flags=
if [ "$(uname -s)" == "Darwin" ] ; then
  cmake_cxx_flags="-I/usr/X11/include -mmacosx-version-min=10.8"
fi

cmake . \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -DVTK_DIR=${VTK_DIR} \
  -DCMAKE_CXX_FLAGS=${cmake_cxx_flags}

make -j8
make install
