#!/usr/bin/env bash
set -ex

#
# the gcc5 mac build output is riddled with innocuous 
#

cmake . -DFS_PACKAGES_DIR="packages" -DBUILD_GUIS=OFF

if [[ "$TRAVIS_OS_NAME" == "osx" ]] ; then
  make 2>&1 | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'
else
  make 
fi
