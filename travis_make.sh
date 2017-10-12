#!/usr/bin/env bash
#set -ex

make -j4 >& build.log
EXITSTATUS=$?
tail -n 20 build.log
exit $EXITSTATUS
