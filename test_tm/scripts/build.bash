#!/bin/bash

# Shell script to build Freesurfer from scratch

# Assumes that the working directory is the 'dev' directory

#   LOGFILE is where to place the output

# Run setup_configure
./setup_configure >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: setup_configure failed"
    exit
fi
echo "NIGHTLY: setup_configure complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE



# Run configure
./configure  --with-boost-dir=/homes/11/rge21/avebury01nfs/boost-1.41.0/ --with-cuda=/usr/local/cuda/3.2.16/cuda/ --enable-fermi-gpu >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: configure failed"
    exit
fi
echo "NIGHTLY: configure complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE



# Do the build
make -j 5 >> $LOGFILE 2>&1
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
