#!/bin/bash

# Unit test for mris_transform. Apply non-linear warp to a surface for
# source and target images with different world matrices.


# Check return status of last command: has to be called immediately
# following the command to test.
function assertpass() {
    if [[ $? -ne 0 ]]; then
        tidyup
        echo FAILED
        exit 1
    fi
}


function tidyup() {
    if [[ ! -z $datadir ]]; then
        rm $datadir/*
        rmdir $datadir
    fi
}


# Set up.
export FREESURFER_HOME=../../distribution
export SUBJECTS_DIR=""
datafile=testdata_01.tar.gz
datadir=$PWD/${datafile%%.*} # Substring before first point.
insurf=$datadir/rh.pial.src
refsurf=$datadir/rh.pial.ref
outsurf=$datadir/rh.pial.out


# Extract test data.
umask 002
tar -zxvf $datafile


# Apply linear warp and check return value of program.
../../mris_transform/mris_transform \
    --src $datadir/src.mgz \
    --dst $datadir/dst.mgz \
    $insurf $datadir/wrp.m3z $outsurf
assertpass


# Compare result against reference surface.
../../mris_diff/mris_diff $outsurf $refsurf
assertpass


tidyup
echo PASSED
exit 0

