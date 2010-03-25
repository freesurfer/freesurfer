#!/bin/bash

# Shell script to run on a cron job

export LOGFILE="${HOME}/nightly.err"

cd ${HOME}/avebury01nfs/dev/

rm $LOGFILE

# Update the directory
cvs update -d . >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: cvs update failed"
    exit
fi

# Do the build
make -j >> $LOGFILE 2>&1

if [ $? -ne 0 ]; then
    echo "NIGHTLY: Build failed"
    exit
fi

# Do the install
make install >> $LOGFILE 2>&1

if [ $? -ne 0 ]; then
    echo "NIGHTLY: Install failed"
    exit
fi

# Run the tests

cd ${HOME}/avebury01nfs/dev/test_tm

tm .

echo "NIGHTLY: Tests complete"