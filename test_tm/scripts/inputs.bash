#!/bin/bash


# Shell script to rsync the gcam/mri etc. inputs to a local drive

# Assumes that
#   SAMPLE_DIR is set to the network share containing the input files
#   SCRATCHDIR is set to the local directory where the files are stored
#   LOGFILE is the file to put place output

mkdir -p $SCRATCHDIR
if [ $? -ne 0 ]; then
    echo "NIGHTLY: mkdir failed"
    exit
fi
echo "NIGHTLY: mkdir complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE

# Get all the input files
rsync -avt --delete $SAMPLE_DIR $SCRATCHDIR  >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: rsync failed"
    exit
fi
echo "NIGHTLY: rsync complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE