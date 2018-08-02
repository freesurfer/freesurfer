#!/usr/bin/env bash

set -ex

#
# script to download the pre-built packages on travis
#

[ "$TRAVIS_OS_NAME" == "linux" ] && tarball="centos7-packages.tar.gz"
[ "$TRAVIS_OS_NAME" == "osx" ]   && tarball="osx10.11-packages.tar.gz"

curl --connect-timeout 8 --retry 60 --retry-delay 3 -O http://surfer.nmr.mgh.harvard.edu/pub/data/fspackages/prebuilt/${tarball}
rm -r packages
tar -xzf ${tarball}
rm ${tarball}
