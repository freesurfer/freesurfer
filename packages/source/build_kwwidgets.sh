#!/usr/bin/env bash

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$1"

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

CMAKE_FLAGS='-DCMAKE_CXX_FLAGS="-I/usr/X11/include -mmacosx-version-min=10.8"'

cmake . -DVTK_DIR=${VTK_DIR} ${CMAKE_FLAGS}

make -j8
make install
