#!/bin/bash

# Shell script to update Freesurfer from CVS

# Assumes that the working directory is the 'dev' directory

#   LOGFILE is where to place the output


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
cvs update -d -P . >> $LOGFILE 2>&1
if [ $? -ne 0 ]; then
    echo "NIGHTLY: cvs update failed"
    exit
fi
echo "NIGHTLY: cvs update complete" >> $LOGFILE

echo >> $LOGFILE
echo >> $LOGFILE
echo >> $LOGFILE