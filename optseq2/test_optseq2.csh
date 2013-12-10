#!/bin/tcsh -f

#############################################################################
#
# Name:    test_optseq2.csh
# Purpose: execute optseq2 and check results against expected results
# Usage:
#
#   test_optseq2.csh
#
#############################################################################

set VERSION='$Id: test_optseq2.csh,v 2.17.4.3 2013/12/10 16:30:47 zkaufman Exp $'

umask 002

set OS=`uname -s`
if ("$OS" == "SunOS") then
  # HACK, dont bother for now
  exit 77
endif

set WD=$PWD

#
# extract testing data
#
gunzip -vc test_data.tar.gz | tar xvf -

#
# extract arch-specific expected results and 
# set location where test data is to be created
#

set EXPECTED=$WD/test_data
cd $EXPECTED
set TEST_DATA=$EXPECTED/emot_testdata.tar.gz

if (! -e $TEST_DATA) then
  echo "Architecture-specific test data file $TEST_DATA does not exist"
  # status 77 means test result is ignored
  exit 77
endif
if (-e emot.sum) rm -Rf emot*
gunzip -vc $TEST_DATA | tar xvf -
if ($status != 0) then
  echo "Failure extracting architecture-specific test data."
  exit 1
endif

set ACTUAL=$WD/tmp
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

set cmd=($WD/test_data/create_optseq2_data.csh $WD/optseq2)
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

#
# remove test data
#
cd $WD
rm -Rf test_data
rm -Rf tmp

echo "test_optseq2 passed all tests"
exit 0
