#!/usr/bin/env bash

set -ex

#
# script to download packages on travis
#

if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    brew cask uninstall oclint
    brew install gcc glib
fi

[ "$TRAVIS_OS_NAME" == "linux" ] && tarball="centos7-packages.tar.gz"
[ "$TRAVIS_OS_NAME" == "osx" ]   && tarball="osx10.11-packages.tar.gz"

cd ..
curl --connect-timeout 8 --retry 60 --retry-delay 3 -O http://surfer.nmr.mgh.harvard.edu/pub/data/fspackages/prebuilt/${tarball}
tar -xzf ${tarball}
rm ${tarball}
