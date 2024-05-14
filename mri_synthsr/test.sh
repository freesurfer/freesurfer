#!/usr/bin/env bash

# trying not to duplicate existing testdata from ../mri_aparc2aseg
# - so do not source common test setup right away

if [ ! -e ../mri_aparc2aseg/testdata.tar.gz ]; then
   echo "Cannot continue w/o existing bert subject data"
   exit 1
fi

# copy in the existing data
if [ -L ../mri_aparc2aseg/testdata.tar.gz ]; then
   # pre-existing bert data in annex is still a soft link to its blob
   slink=`ls -l ../mri_aparc2aseg/testdata.tar.gz | sed 's;^.*-> ;;'`
   ls -l $slink
   rm -f testdata.tar.gz
   cp -p -f $slink testdata.tar.gz
   ls -l testdata.tar.gz
elif [ -f ../mri_aparc2aseg/testdata.tar.gz ]; then
   # already unlocked from the annex
   rm -f testdata.tar.gz
   cp -p -f ../mri_aparc2aseg/testdata.tar.gz .
else
   echo "Cannot determine status of bert subject data testdata.tar.gz"
   exit 1
fi

# use_cmake_test_framework="false"
use_cmake_test_framework="true"

if [ "$use_cmake_test_framework" == "true" ]; then
   source "$(dirname $0)/../test.sh"
   # running with 4 threads should not put too much stress on CPU/memory usage
   # - given other tests are running (in parallel) during nightly builds
   # which mri_synthsr
   test_command mri_synthsr --i $SUBJECTS_DIR/bert/mri --o bert.out --threads 4
   rm -f testdata.tar.gz
elif [ "$use_cmake_test_framework" == "false" ]; then
   # run test standalone, e.g., debugging
   rm -f testdata.tar
   gunzip testdata.tar.gz
   rm -rf testdata
   tar xpf testdata.tar
   export FREESURFER_HOME=/usr/local/freesurfer/dev
   source $FREESURFER_HOME/SetUpFreeSurfer.sh
   # (cd testdata && date && which mri_synthsr)
   (cd testdata && date && bash -x `which mri_synthsr` --i $SUBJECTS_DIR/bert/mri --o bert.out --threads 4 && date)
else
   echo "Nothing selected to run."
   exit 1
fi

