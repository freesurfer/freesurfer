#!/usr/bin/env bash

# Even on OS X we should build with the "linux-g++" configuration. The only real
# difference is that the linux config uses ar instead of libtool, and the libtool flags
# (specifically -static) cause problems since it's not expecting glibtool (which we usually
# use on mac). Another note: the only differences between ANN 1.1.1 (which we use at Martinos) 
# and 1.1.2 are a couple very small changes to hide modern gcc errors

set -e

if [ "$#" != "1" ] ; then echo "error: usage: build.sh <prefix>" && exit 1 ; fi
INSTALL_DIR="$(realpath $1)"

cd ann_1.1.2

make linux-g++
rm -rf ${INSTALL_DIR}/lib ${INSTALL_DIR}/include
mv lib include ${INSTALL_DIR}
