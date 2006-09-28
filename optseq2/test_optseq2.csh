#!/bin/tcsh -ef

#############################################################################
#
# Name:    test_optseq2.csh
# Purpose: execute optseq2 and check results against expected results
# Usage:
#
#   test_optseq2.csh
#
#############################################################################

set VERSION='$Id: test_optseq2.csh,v 2.9 2006/09/28 06:34:05 nicks Exp $'

umask 002

set WD=$PWD

#
# gunzip arch-specific expected results and 
# set location where test data is to be created
#

set EXPECTED=$WD/test_data
cd $EXPECTED
set PROC=`uname -p`
set TEST_DATA=$EXPECTED/$PROC-emot.tar.gz
if (! -e $TEST_DATA) then
  echo "Architecture-specific test data file $TEST_DATA does not exist"
  # status 77 means test result is ignored
  exit 77
endif
if (-e emot.sum) rm -Rf emot*
set cmd=(tar zxvf $TEST_DATA)
echo $cmd
$cmd
if ($status != 0) then
  echo "Failure extracting architecture-specific test data."
  exit 1
endif
chmod 666 emot*

set ACTUAL=/tmp/optseq2
if (-e $ACTUAL) rm -Rf $ACTUAL
mkdir -p $ACTUAL
if (! -e $ACTUAL) then
    echo "test_optseq2 FAILED to create directory $ACTUAL"
    exit 1
endif
chgrp -R fsdev $ACTUAL

#
# run optseq2 using typical input args (and the same as those used
# to produce the expected results)
#

set cmd=($WD/test_data/create_optseq2_data.csh $WD/optseq2)
cd $ACTUAL
echo $cmd
$cmd
if ($status != 0) then
  echo "create_optseq2_data FAILED"
  chmod -R 777 $ACTUAL
  exit 1
endif
chmod -R 777 $ACTUAL

#
# compare expected results with actual (produced) results
#

set TEST_FILES=( emot-001.par emot-002.par emot-003.par emot-004.par emot.sum )
set filter=($WD/test_data/filter_optseq2_data.csh)
foreach tstfile ($TEST_FILES)
  # first remove irrelevant lines of data that may change from run-to-run
  set cmd=($filter $EXPECTED/$tstfile)
  $cmd
  set cmd=($filter $ACTUAL/$tstfile)
  $cmd
  # now compare expected with observed data
  set cmd=(diff $EXPECTED/$tstfile $ACTUAL/$tstfile);
  echo $cmd
  $cmd
  set diff_status=$status
  if ($diff_status != 0) then
    echo "$cmd FAILED (exit status=$diff_status)"
    exit 1
  endif
end

# cleanup for the next run
chgrp -R fsdev $EXPECTED
rm -Rf $ACTUAL

echo "test_optseq2 passed all tests"
exit 0
