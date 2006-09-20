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

set VERSION='$Id: test_optseq2.csh,v 2.1 2006/09/20 19:53:39 nicks Exp $'

umask 002

#
# point to expected results and location where test data is to be created
#

set EXPECTED=$PWD/test_data
set ACTUAL=/tmp/optseq2
if (-e $ACTUAL) rm -Rf $ACTUAL
mkdir -p $ACTUAL
if (! -e $ACTUAL) then
    echo "test_optseq2 FAILED to create directory $ACTUAL"
    exit 1
endif

#
# run optseq2 using typical input args (and the same as those used
# to produce the expected results)
#

set cmd=($PWD/test_data/create_optseq2_data.csh)
cd $ACTUAL
echo $cmd
$cmd
if ($status != 0) then
  echo "create_optseq2_data FAILED"
  exit 1
endif

#
# compare expected results with actual (produced) results
#

set TEST_FILES=( emot-001.par emot-002.par emot-003.par emot-004.par )
foreach tstfile ($TEST_FILES)
  set cmd=(diff $EXPECTED/$tstfile $ACTUAL/$tstfile);
  echo $cmd
  $cmd
  set diff_status=$status
  if ($diff_status != 0) then
    echo "$cmd FAILED (exit status=$diff_status)"
    exit 1
  endif
end

if (-e $ACTUAL) rm -Rf $ACTUAL

exit 0
