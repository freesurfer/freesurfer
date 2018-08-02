#!/usr/bin/env bash

set -ex
set -o pipefail

#
# the gcc5 mac output is riddled with innocuous deprecation warnings -
# we'll have to deal with these until gcc6, but for now, we should filter the output
# so the travis log doesn't get filled up
#

cmake . -DFS_PACKAGES_DIR="packages" -DBUILD_GUIS=OFF

if [[ "$TRAVIS_OS_NAME" == "osx" ]] ; then
  make 2>&1 | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'
else
  make 
fi
