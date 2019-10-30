#!/usr/bin/env bash

set -ex

#
# script to configure and build on travis
#

cmake . -DFS_PACKAGES_DIR="$(pwd)/../packages" -DBUILD_GUIS=OFF
make -j4
