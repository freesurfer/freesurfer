#!/usr/bin/env bash

set -ex
set -o pipefail

#
# the gcc5 mac output is riddled with innocuous deprecation warnings -
# we'll have to deal with these until gcc6, but for now, we should filter the output
# so the travis log doesn't get filled up
#

nm -D /usr/lib/x86_64-linux-gnu/libGLU.so | grep gluErrorString
nm -D /usr/lib/x86_64-linux-gnu/libGLU.so | grep gluBuild2DMipmaps
nm -D /usr/lib/x86_64-linux-gnu/libGL.so  | grep gluErrorString
nm -D /usr/lib/x86_64-linux-gnu/libGL.so  | grep gluBuild2DMipmaps

cmake . -DFS_PACKAGES_DIR="packages" -DBUILD_GUIS=OFF

if [[ "$TRAVIS_OS_NAME" == "osx" ]] ; then
  make 2>&1 | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'
else
  make mris2rgb VERBOSE=1
fi
