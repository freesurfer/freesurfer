#!/usr/bin/env bash
set -ex

#
# script to download the pre-built packages on travis
#

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  curl --verbose --connect-timeout 8 --retry 60 --retry-delay 3 -O http://surfer.nmr.mgh.harvard.edu/pub/dist/fs_supportlibs/prebuilt/centos6_x86_64/centos6-x86_64-packages.tar.gz
  tar -xzf centos6-x86_64-packages.tar.gz 
  rm centos6-x86_64-packages.tar.gz
  mv centos6-x86_64-packages build-packages
  cd build-packages
  ./setup.sh
  cd ..
elif [[ "$TRAVIS_OS_NAME" == "osx" ]];   then
  mkdir build-packages
  cd build-packages
  curl --verbose --connect-timeout 8 --retry 60 --retry-delay 3 -O http://surfer.nmr.mgh.harvard.edu/pub/dist/fs_supportlibs/prebuilt/OSX/osx-lion-packages.tar.gz
  tar -xzf osx-lion-packages.tar.gz
  rm osx-lion-packages.tar.gz
  cd ..
fi
