#!/bin/bash

# Shell script to run on a cron job

# Where to keep a log
export LOGFILE="${HOME}/nightly.err"

# Location of the nightly build directory
export NIGHTLYDIR="${HOME}/avebury01nfs/dev/"

# 'touch' is to prevent complaints if file doesn't exist
touch $LOGFILE
rm $LOGFILE

# Make sure environment is set up
export PS1="--"
source $HOME/.bashrc >> $LOGFILE 2>&1

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE

# ---------------------------

# Go to the build directory
cd $NIGHTLYDIR
pwd

# Clean it
make clean >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: make clean failed"
    exit
fi
echo "NIGHTLY: make clean complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE

# Update the directory
cvs update -d . >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: cvs update failed"
    exit
fi
echo "NIGHTLY: cvs update complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE

# Do the build
make -j >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: Build failed"
    exit
fi
echo "NIGHTLY: Build complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE

# Do the install
make install >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: Install failed"
    exit
fi
echo "NIGHTLY: Installation complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE


# Run the tests

cd $NIGHTLYDIR/test_tm
echo
pwd
echo

echo "Running tm"
echo
tm .
echo

echo "NIGHTLY: Tests complete"