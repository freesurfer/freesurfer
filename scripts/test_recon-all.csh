#! /bin/tcsh -f

#
# Name:    test_recon-all.csh
# Purpose: runs recon-all on a reference subject, then checks results
# Usage:
#
#   test_recon-all.csh [-rsd <reference subj source dir>] \
#                      [-rs <reference subjid>] \
#                      [-tsd <test subject dest dir>] \
#                      [-ts <test subjid>] \
#                      [-fshome <FREESURFER_HOME>]
#
#   the defaults are:
#     <reference subj source dir> = 
# /space/freesurfer/subjects/test/weekly_test/subjects/`uname -p`
#     <reference subjid> = bert
#     <test subject dest dir> = /tmp
#     <test subjid> = bert
#     <FREESURFER_HOME> = /usr/local/freesurfer/stable
#
#   the utilities run by this script include:
#     recon-all
#     mri_diff
#     mri_compute_seg_overlap
#     mris_diff
#     mri_surf2surf
#     mris_compute_parc_overlap
#     diff
#     asegstatsdiff
#     aparcstatsdiff
#
# Original Author: Nick Schmansky
#
# Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
#
# Terms and conditions for use, reproduction, distribution and contribution
# are found in the 'FreeSurfer Software License Agreement' contained
# in the file 'LICENSE' found in the FreeSurfer distribution, and here:
#
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
#
# Reporting: freesurfer@nmr.mgh.harvard.edu
#
#


set VERSION='test_recon-all.csh @FS_VERSION@'

set MAIL_LIST=(zkaufman@nmr.mgh.harvard.edu)
# failure mailing list:
set FMAIL_LIST=(zkaufman@nmr.mgh.harvard.edu)

limit coredumpsize unlimited

umask 002

# for debug, setenv DONT_RUNIT to disable command execution
set RunIt=1
if ($?DONT_RUNIT) set RunIt=0
# for debug, setenv SET_ECHO_1 to show all command execution
unsetenv echo
if ($?SET_ECHO_1) set echo=1

#
# setup defaults:
set PROC=`uname -p`
# reference subject: per machine processor type
set SUBJ_REF_DIR=/space/freesurfer/subjects/test/weekly_test/subjects/$PROC
# name of reference subject:
set REF_SUBJ=bert
# this is where the recon-all results will go (the test subject):
setenv SUBJECTS_DIR /tmp
# name to give to test subject:
set TEST_SUBJ=bert
# freesurfer enviro to use
setenv FREESURFER_HOME /usr/local/freesurfer/stable

#
# override defaults via command-line args:
set cmdline = ($argv);
while( $#argv != 0 )

  set flag = $argv[1]; shift;

  switch($flag)

    case "-rsd":
      if ( $#argv < 1) goto arg1err;
      set SUBJ_REF_DIR=$argv[1]; shift;
      breaksw

    case "-rs":
      if ( $#argv < 1) goto arg1err;
      set REF_SUBJ=$argv[1]; shift;
      set REF_SUBJ=`basename $REF_SUBJ`; # removes trailing /
      breaksw

    case "-tsd":
      if ( $#argv < 1) goto arg1err;
      setenv SUBJECTS_DIR $argv[1]; shift;
      breaksw

    case "-ts":
      if ( $#argv < 1) goto arg1err;
      set TEST_SUBJ=$argv[1]; shift;
      set TEST_SUBJ=`basename $TEST_SUBJ`; # removes trailing /
      breaksw

    case "-fshome":
      if ( $#argv < 1) goto arg1err;
      setenv FREESURFER_HOME $argv[1]; shift;
      breaksw

    case "-parallel":
      setenv PARALLEL 1;
      breaksw

    case "-norecon":
      setenv SKIP_RECON 1
      breaksw

    case "-help":
    case "--help":
      echo "test_recon-all.csh [-rsd <reference subj source dir>]"
      echo "                   [-rs <reference subjid>]"
      echo "                   [-tsd <test subject dest dir>]"
      echo "                   [-ts <test subjid>]"
      echo "                   [-fshome <FREESURFER_HOME>]"
      echo "                   [-norecon]"
      echo "\nthe defaults are:"
      echo "  <reference subj source dir> = /space/freesurfer/subjects/test/weekly_test/subjects/$PROC"
      echo "  <reference subjid> = bert"
      echo "  <test subject dest dir> = /tmp"
      echo "  <test subjid> = bert"
      echo "  <FREESURFER_HOME> = /usr/local/freesurfer/stable"
      echo "\nthe utilities run by this script include:"
      echo "  recon-all (unless -norecon is specified)"
      echo "  mri_diff"
      echo "  mri_compute_seg_overlap"
      echo "  mris_diff"
      echo "  mri_surf2surf"
      echo "  mris_compute_parc_overlap"
      echo "  diff"
      echo "  asegstatsdiff"
      echo "  aparcstatsdiff"
      exit 1
      breaksw

    default:
      echo ERROR: Flag $flag unrecognized.
      echo $cmdline
      exit 1
      breaksw
  endsw
end

echo "SUBJ_REF_DIR:    $SUBJ_REF_DIR"
echo "REF_SUBJ:        $REF_SUBJ"
echo "SUBJECTS_DIR:    $SUBJECTS_DIR"
echo "TEST_SUBJ:       $TEST_SUBJ"
echo "FREESURFER_HOME: $FREESURFER_HOME"

#
# first check if the prior test run failed, in which case don't bother
# running another test run until somebody figures out why it failed.
# the flag file 'test_recon-all_FAILED' will  need to be deleted manually.
if (-e $SUBJECTS_DIR/test_recon-all_FAILED) then
    set msg="Prior $PROC test_recon-all FAILED on $HOST! Fix it!"
    echo ${msg}
    mail -v -s "${msg}" $MAIL_LIST < /dev/null > /dev/null
    exit 1
endif

#
# next check if recon-all is already running, as indicated by
# existence of the subject dir, which is moved upon successful
# completion of recon-all, and so should not exist
if (-e $SUBJECTS_DIR/$TEST_SUBJ) then
  if ($?SKIP_RECON) goto continue1
    set msg="$PROC test_recon-all already running on $HOST !"
    echo ${msg}
    mail -v -s "${msg}" $MAIL_LIST < /dev/null > /dev/null
    exit 1
  continue1:
endif

#
# and create a symlink in our SUBJECTS_DIR naming the reference subject
# as 'ref_subj', ensuring that both the reference subject and the test
# subject (as yet created) are in the same SUBJECTS_DIR, as some utilities
# assume that all subjects are in the same SUBJECTS_DIR
rm -f $SUBJECTS_DIR/ref_subj
ln -s $SUBJ_REF_DIR/$REF_SUBJ $SUBJECTS_DIR/ref_subj
if ($status) then
    set msg="test_recon-all FAILED to create ref_subj symlink!"
    echo ${msg}
    mail -v -s "${msg}" $MAIL_LIST < /dev/null > /dev/null
    exit 1
endif


#
# setup and begin logging
#
set LOG_DIR=/space/freesurfer/subjects/test/weekly_test/logs/$PROC/$TEST_SUBJ
mkdir -p $LOG_DIR
set OUTPUTF=$LOG_DIR/test_recon-all.txt
set RECON_LOG=$SUBJECTS_DIR/$TEST_SUBJ/scripts/recon-all.log
echo "$PROC recon-all test for $TEST_SUBJ" >& $OUTPUTF
set BEGIN_TIME=`date`
echo "Version: $VERSION" >>& $OUTPUTF
if ($#argv) echo "args: $argv" >>& $OUTPUTF
echo "Start: $BEGIN_TIME" >>& $OUTPUTF

source $FREESURFER_HOME/SetUpFreeSurfer.csh

# setenv SKIP_RECON to skip recon-all and procede to compare
# previously completed results against the reference data:
if ($?SKIP_RECON) goto compare
echo "Starting recon-all in five seconds... (Ctrl-C to stop)"
sleep 5

# check if mni tools were built against /usr/pubsw/bin/perl 
# instead of /usr/bin/perl
grep pubsw `which nu_correct`
if ( ! $status ) then
    echo "***FAILED :: MNI tools built using /usr/pubsw/perl!"  >>& $OUTPUTF
    mail -v -s "test_recon-all -all FAILED: MNI tools build using /usr/pubsw/perl!" $FMAIL_LIST < $RECON_LOG
    cp -f $RECON_LOG $LOG_DIR/
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    exit 1
endif

#
# run recon-all
#

# setup the libsafe buffer overflow and format string violation detector
if (-e /lib/libsafe.so.2) setenv LD_PRELOAD /lib/libsafe.so.2

#
# gather-up possible input volumes
set INVOL_LIST=(001.mgz 002.mgz 003.mgz 004.mgz 005.mgz)
set INVOL_LIST=($INVOL_LIST 006.mgz 007.mgz 008.mgz 009.mgz 010.mgz)
set INVOL_LIST=($INVOL_LIST 011.mgz 012.mgz 013.mgz 014.mgz 015.mgz)
set INVOL_LIST=($INVOL_LIST 016.mgz 017.mgz 018.mgz 019.mgz 020.mgz)
set INVOL_LIST=($INVOL_LIST 021.mgz 022.mgz 023.mgz 024.mgz 025.mgz)
set INVOL_LIST=($INVOL_LIST 026.mgz 027.mgz 028.mgz 029.mgz 030.mgz)
set INVOL_LIST=($INVOL_LIST 031.mgz 032.mgz 033.mgz 034.mgz 035.mgz)
set INVOL=()
set echo=1
foreach invol ($INVOL_LIST)
  set involfull=$SUBJECTS_DIR/ref_subj/mri/orig/$invol
  if (-e $involfull) set INVOL=($INVOL -i $involfull)
end
if ($#INVOL == "0") then
    echo "***FAILED :: no input volumes found"  >>& $OUTPUTF
    mail -v -s "test_recon-all -all FAILED: no input volumes found" $FMAIL_LIST < $RECON_LOG
    cp -f $RECON_LOG $LOG_DIR/
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    exit 1
endif

# run a copy of recon-all, since the nightly build can clobber recon-all
# if it changed
cp -f `which recon-all` /tmp
if ($status) exit 1

# set recon-all command and run...
set cmd=(/tmp/recon-all)
set cmd=($cmd -s $TEST_SUBJ $INVOL)
set cmd=($cmd -all -debug -clean -norandomness -allowcoredump -time);
if ($?PARALLEL) then
  set cmd=($cmd -parallel)
endif
echo $cmd
if ($RunIt) then
  cd $SUBJECTS_DIR
  # recon-all: this will take some 40 hours to run...
  /usr/bin/time $cmd >& $SUBJECTS_DIR/recon-all.log.txt
  if ($status != 0) then
    echo "***FAILED :: $PROC recon-all -all"  >>& $OUTPUTF
    mail -v -s "test_recon-all -all FAILED on $PROC" $FMAIL_LIST < $RECON_LOG
    cp -f $RECON_LOG $LOG_DIR/
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/recon-all.log.txt
    exit 1
  else
    set CURRENT_TIME=`date`
    set ELAPSED=`grep elapsed $SUBJECTS_DIR/recon-all.log.txt | grep CPU | awk '{print $3}'`
    echo "   pass :: recon-all -all (Finish: $CURRENT_TIME, ${ELAPSED})" >>& $OUTPUTF
    cp -f $RECON_LOG $LOG_DIR/
    chmod a+w $SUBJECTS_DIR/recon-all.log.txt
  endif
endif # ($RunIt)

#
# check if recon-all exited without failure code, but did not actually finish
# 
if (-e $SUBJECTS_DIR/$TEST_SUBJ/scripts/IsRunning.lh+rh) then
    echo "$PROC test_recon-all FAILED"  >>& $OUTPUTF
    echo "IsRunning.lh+rh flag exists: recon-all did not finish" >>& $OUTPUTF
    echo "Check recon-all.log to pinpoint failure" >>& $OUTPUTF
    mail -v -s "test_recon-all FAILED on $PROC" $FMAIL_LIST < $OUTPUTF
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    exit 1
endif


#
# compare resulting volumes and surfaces with reference data (prior recon)
#

compare:

if (! -e $SUBJECTS_DIR/$TEST_SUBJ) then
    echo "$PROC test_recon-all FAILED"  >>& $OUTPUTF
    echo "missing $SUBJECTS_DIR/$TEST_SUBJ" >>& $OUTPUTF
    mail -v -s "test_recon-all FAILED on $PROC" $FMAIL_LIST < $OUTPUTF
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    exit 1
endif

set TEST_VOLUMES=(rawavg.mgz orig.mgz nu.mgz T1.mgz brainmask.mgz \
norm.mgz aseg.mgz brain.mgz wm.mgz filled.mgz aparc+aseg.mgz \
lh.ribbon.mgz rh.ribbon.mgz)

foreach tstvol ($TEST_VOLUMES)
  set REF_VOL  = $SUBJECTS_DIR/ref_subj/mri/$tstvol
  set TST_VOL  = $SUBJECTS_DIR/$TEST_SUBJ/mri/$tstvol
  set DIFF_VOL = $LOG_DIR/mri_diff-$tstvol
  set MRIDIFFF = $LOG_DIR/mri_diff-$tstvol.txt
  set cmd=(mri_diff --debug --thresh 0 --log $MRIDIFFF \
           $REF_VOL $TST_VOL --diff $DIFF_VOL);
  echo $cmd
  if ($RunIt) then
    if (-e $MRIDIFFF) rm -f $MRIDIFFF
    $cmd
    set mri_diff_status=$status
    if ($mri_diff_status != 0) then
      printf "***FAILED :: mri_diff $tstvol (exit status=$mri_diff_status)\n" \
           >>& $OUTPUTF
      setenv FOUND_ERROR 1
      # continue running tests
    else
      printf "   pass :: mri_diff $tstvol\n" >>& $OUTPUTF
    endif
    if (-e $MRIDIFFF) chmod g+rw $MRIDIFFF
    if (-e $DIFF_VOL) chmod g+rw $DIFF_VOL

    # squeeze-in a check of the transform matrices, where appropriate
    if ("$tstvol" == "brainmask.mgz") then
      foreach xform (talairach.xfm talairach.lta talairach_with_skull.lta)
        # filter-out time-stamp and filename from .lta files:
        grep -v created $SUBJECTS_DIR/ref_subj/mri/transforms/$xform \
            > $SUBJECTS_DIR/ref.tmp.$xform
        grep -v filename $SUBJECTS_DIR/ref.tmp.$xform \
            > $SUBJECTS_DIR/ref.$xform
        rm $SUBJECTS_DIR/ref.tmp.$xform
        grep -v created $SUBJECTS_DIR/$TEST_SUBJ/mri/transforms/$xform \
            > $SUBJECTS_DIR/tst.tmp.$xform
        grep -v filename $SUBJECTS_DIR/tst.tmp.$xform \
            > $SUBJECTS_DIR/tst.$xform
        rm $SUBJECTS_DIR/tst.tmp.$xform
        set cmd=(diff $SUBJECTS_DIR/ref.$xform $SUBJECTS_DIR/tst.$xform)
        $cmd >& $LOG_DIR/diff-$xform.txt
        set diff_status=$status
        if ($diff_status != 0) then
          printf "***FAILED :: diff $xform\n" >>& $OUTPUTF
          setenv FOUND_ERROR 1
          # continue running tests
        else
          printf "   pass :: diff $xform\n" >>& $OUTPUTF
        endif
      end
    endif

  endif # ($RunIt)
end

set REF_VOL = $SUBJECTS_DIR/ref_subj/mri/aseg.mgz
set TST_VOL = $SUBJECTS_DIR/$TEST_SUBJ/mri/aseg.mgz
set cmd=(mri_compute_seg_overlap $REF_VOL $TST_VOL)
echo $cmd
if ($RunIt) then
  $cmd >& $LOG_DIR/mri_compute_seg_overlap.txt
  grep "Overall subcortical dice = 1.0000" $LOG_DIR/mri_compute_seg_overlap.txt
  if ($status) then
    printf "***FAILED :: mri_compute_seg_overlap aseg.mgz (Dice != 1)\n" \
        >>& $OUTPUTF
    setenv FOUND_ERROR 1
    # continue running tests
  else
    printf "   pass :: mri_compute_seg_overlap aseg.mgz\n" >>& $OUTPUTF
  endif
endif

set TEST_SURFACES=(orig.nofix smoothwm.nofix inflated.nofix qsphere.nofix \
orig smoothwm inflated white pial sphere sphere.reg)
# note: mris_diff cannot check jacobian_white)
set TEST_HEMIS=(rh lh)

foreach hemi ($TEST_HEMIS)
  foreach tstsurf ($TEST_SURFACES)
    set REF_SURF = $SUBJECTS_DIR/ref_subj/surf/$hemi.$tstsurf
    set TST_SURF = $SUBJECTS_DIR/$TEST_SUBJ/surf/$hemi.$tstsurf
    set MRISDIFFF=$LOG_DIR/mris_diff-$hemi.$tstsurf.txt
    set cmd=(mris_diff --debug --thresh 0 --maxerrs 1000)
    set cmd=($cmd $REF_SURF $TST_SURF)
    echo $cmd
    if ($RunIt) then
      if (-e $MRISDIFFF) rm -f $MRISDIFFF
      $cmd >& $MRISDIFFF
      set mris_diff_status=$status
      if ($mris_diff_status != 0) then
        printf "***FAILED :: mris_diff $hemi.$tstsurf (exit status=$mris_diff_status)\n"  >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # continue running tests
      else
        printf "   pass :: mris_diff $hemi.$tstsurf\n" >>& $OUTPUTF
      endif
      if (-e $MRISDIFFF) chmod g+rw $MRISDIFFF
    endif # ($RunIt)
  end
end


set TEST_LABELS=(cortex.label)
set TEST_HEMIS=(rh lh)

foreach hemi ($TEST_HEMIS)
  foreach labelname ($TEST_LABELS)
    set DIFFF=$LOG_DIR/diff-$hemi.$labelname.txt
    set cmd=(diff)
    set cmd=($cmd $SUBJECTS_DIR/ref_subj/label/${hemi}.${labelname})
    set cmd=($cmd $SUBJECTS_DIR/${TEST_SUBJ}/label/${hemi}.${labelname})
    echo $cmd
    if ($RunIt) then
      if (-e $DIFFF) rm -f $DIFFF
      $cmd >& $DIFFF
      set diff_status=$status
      if ($diff_status != 0) then
        printf "***FAILED :: diff $hemi.$labelname (exit status=$diff_status)\n"  >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # continue running tests
      else
        printf "   pass :: diff $hemi.$labelname\n" >>& $OUTPUTF
      endif
      if (-e $DIFFF) chmod g+rw $DIFFF
    endif # ($RunIt)
  end
end


set TEST_CURVS=(curv curv.pial sulc thickness area area.pial volume)
set TEST_HEMIS=(rh lh)

foreach hemi ($TEST_HEMIS)
  foreach curvname ($TEST_CURVS)
    set MRISDIFFF=$LOG_DIR/mris_diff-$hemi.$curvname.txt
    set cmd=(mris_diff --debug --thresh 0 --maxerrs 1000)
    set cmd=($cmd --sd1 $SUBJECTS_DIR --s1 ref_subj)
    set cmd=($cmd --sd2 $SUBJECTS_DIR --s2 $TEST_SUBJ)
    set cmd=($cmd --hemi $hemi)
    set cmd=($cmd --curv $curvname)
    echo $cmd
    if ($RunIt) then
      if (-e $MRISDIFFF) rm -f $MRISDIFFF
      $cmd >& $MRISDIFFF
      set mris_diff_status=$status
      if ($mris_diff_status != 0) then
        printf "***FAILED :: mris_diff $hemi.$curvname (exit status=$mris_diff_status)\n"  >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # create a difference file that can be loaded as an overlay in
        # tksurfer, showing exactly where the differences are found
        set REF_CURV = $SUBJECTS_DIR/ref_subj/surf/$hemi.$curvname
        set TST_CURV = $SUBJECTS_DIR/$TEST_SUBJ/surf/$hemi.$curvname
        set cmd=(mri_diff --diff $LOG_DIR/mri_diff-$hemi.$curvname.mgz)
        set cmd=($cmd $REF_CURV $TST_CURV)
        echo $cmd
        $cmd
        # continue running tests
      else
        printf "   pass :: mris_diff $hemi.$curvname\n" >>& $OUTPUTF
      endif
      if (-e $MRISDIFFF) chmod g+rw $MRISDIFFF
    endif # ($RunIt)
  end
end

compare2:

set TEST_APARCS=(aparc.a2009s aparc)
set TEST_HEMIS=(rh lh)

foreach hemi ($TEST_HEMIS)
  foreach aparcname ($TEST_APARCS)
    set MRISDIFFF=$LOG_DIR/mris_diff-$hemi.$aparcname.txt
    set cmd=(mris_diff --debug --maxerrs 1000)
    set cmd=($cmd --sd1 $SUBJECTS_DIR --s1 ref_subj)
    set cmd=($cmd --sd2 $SUBJECTS_DIR --s2 $TEST_SUBJ)
    set cmd=($cmd --hemi $hemi)
    set cmd=($cmd --aparc $aparcname)
    echo $cmd
    if ($RunIt) then
      if (-e $MRISDIFFF) rm -f $MRISDIFFF
      $cmd >& $MRISDIFFF
      set mris_diff_status=$status
      if ($mris_diff_status != 0) then
        printf "***FAILED :: mris_diff $hemi.$aparcname (exit status=$mris_diff_status)\n"  >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # continue running tests
      else
        printf "   pass :: mris_diff $hemi.$aparcname\n" >>& $OUTPUTF
      endif
      if (-e $MRISDIFFF) chmod g+rw $MRISDIFFF
    endif # if ($RunIt)

#
# to get a handle on changes in parcellation, first run mri_surf2surf
# to map the reference annotations onto the test surface, then run
# mris_compute_parc_overlap to calculate a Dice coefficent, and more
# usefuly, create some .annot files which show which areas are different
    set cmd=(mri_surf2surf)
    set cmd=($cmd --srcsubject ref_subj)
    set cmd=($cmd --trgsubject $TEST_SUBJ)
    set cmd=($cmd --hemi $hemi)
    set cmd=($cmd --sval-annot \
        $SUBJECTS_DIR/ref_subj/label/$hemi.$aparcname.annot)
    set cmd=($cmd --tval \
        $SUBJECTS_DIR/$TEST_SUBJ/label/$hemi.${REF_SUBJ}_ref.$aparcname.annot)
    echo $cmd
    if ($RunIt) then
      set S2SF=$LOG_DIR/mri_surf2surf-$hemi.$aparcname.txt
      if (-e $S2SF) rm -f $S2SF
      $cmd >& $S2SF
      if ($status) then
        printf "***FAILED :: mri_surf2surf $hemi.$aparcname\n" >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # continue running tests
      else
        printf "   pass :: mri_surf2surf $hemi.$aparcname\n" >>& $OUTPUTF
      endif
      if (-e $S2SF) chmod g+rw $S2SF
    endif # if ($RunIt)

    set cmd=(mris_compute_parc_overlap)
    set cmd=($cmd --sd $SUBJECTS_DIR)
    set cmd=($cmd --s $TEST_SUBJ)
    set cmd=($cmd --hemi $hemi)
    set cmd=($cmd --annot1 $aparcname)
    set cmd=($cmd --annot2 ${REF_SUBJ}_ref.$aparcname)
    set cmd=($cmd --debug-overlap)
    echo $cmd
    if ($RunIt) then
      set MCPOF=$LOG_DIR/mris_comp_parc_olap-$hemi.$aparcname.txt
      if (-e $MCPOF) rm -f $MCPOF
      $cmd >& $MCPOF
      if ($status) then
        printf "***FAILED :: mris_compute_parc_overlap $hemi.$aparcname\n" \
            >>& $OUTPUTF
        setenv FOUND_ERROR 1
        # continue running tests
      else
        printf "   pass :: mris_compute_parc_overlap $hemi.$aparcname\n" \
            >>& $OUTPUTF
      endif
      if (-e $MCPOF) chmod g+rw $MCPOF
    endif # if ($RunIt)

  end # foreach aparcname ($TEST_APARCS)
end # foreach hemi ($TEST_HEMIS)


set STATS_FILES=(aseg.stats \
lh.aparc.a2009s.stats lh.aparc.stats \
rh.aparc.a2009s.stats rh.aparc.stats wmparc.stats)

foreach statfile ($STATS_FILES)
  set REF_STAT   = $SUBJECTS_DIR/ref_subj/stats/$statfile
  set TST_STAT   = $SUBJECTS_DIR/$TEST_SUBJ/stats/$statfile
  set STATSDIFFF = $LOG_DIR/diff-$statfile.txt

  # grep-out stuff that legitimately changes from run-to-run
  set cmd1a=(grep -v "TimeStamp" $REF_STAT)
  set cmd1b=(grep -v "CreationTime" $SUBJECTS_DIR/ref.stats)
  set cmd1c=(grep -v "cvs_version" $SUBJECTS_DIR/ref.stats)
  set cmd1d=(grep -v "${REF_SUBJ}" $SUBJECTS_DIR/ref.stats)
  set cmd1e=(grep -v "hostname" $SUBJECTS_DIR/ref.stats)
  set cmd1f=(grep -v "SUBJECTS_DIR" $SUBJECTS_DIR/ref.stats)
  set cmd1g=(grep -v "machine" $SUBJECTS_DIR/ref.stats)
  set cmd1h=(grep -v "ColorTable" $SUBJECTS_DIR/ref.stats)
  echo $cmd1a
  echo $cmd1b
  echo $cmd1c
  echo $cmd1d
  echo $cmd1e
  echo $cmd1f
  echo $cmd1g
  echo $cmd1h

  set cmd2a=(grep -v "TimeStamp" $TST_STAT)
  set cmd2b=(grep -v "CreationTime" $SUBJECTS_DIR/tst.stats)
  set cmd2c=(grep -v "cvs_version" $SUBJECTS_DIR/tst.stats)
  set cmd2d=(grep -v "${TEST_SUBJ}" $SUBJECTS_DIR/tst.stats)
  set cmd2e=(grep -v "hostname" $SUBJECTS_DIR/tst.stats)
  set cmd2f=(grep -v "SUBJECTS_DIR" $SUBJECTS_DIR/tst.stats)
  set cmd2g=(grep -v "machine" $SUBJECTS_DIR/ref.stats)
  set cmd2h=(grep -v "ColorTable" $SUBJECTS_DIR/ref.stats)
  echo $cmd2a
  echo $cmd2b
  echo $cmd2c
  echo $cmd2d
  echo $cmd2e
  echo $cmd2f
  echo $cmd2g
  echo $cmd2h

  set cmd3=(diff $SUBJECTS_DIR/ref.$statfile $SUBJECTS_DIR/tst.$statfile)
  echo $cmd3

  if ($RunIt) then
    if (-e $STATSDIFFF) rm -f $STATSDIFFF
    $cmd1a >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.stats
    $cmd1b >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.stats
    $cmd1c >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.stats 
    $cmd1d >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.stats
    $cmd1e >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.stats
    $cmd1f >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.$statfile
    $cmd1g >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.$statfile
    $cmd1h >& $SUBJECTS_DIR/ref.stats.tmp
    mv $SUBJECTS_DIR/ref.stats.tmp $SUBJECTS_DIR/ref.$statfile

    $cmd2a >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.stats
    $cmd2b >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.stats
    $cmd2c >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.stats
    $cmd2d >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.stats
    $cmd2e >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.stats
    $cmd2f >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.$statfile
    $cmd2g >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.$statfile
    $cmd2h >& $SUBJECTS_DIR/tst.stats.tmp
    mv $SUBJECTS_DIR/tst.stats.tmp $SUBJECTS_DIR/tst.$statfile

    $cmd3 >& $STATSDIFFF
    set stats_diff_status=$status
    if ($stats_diff_status != 0) then
      printf "***FAILED :: diff $statfile (exit status=$stats_diff_status)\n" \
           >>& $OUTPUTF
      setenv FOUND_ERROR 1
      # continue running tests
    else
      printf "   pass :: diff $statfile\n" >>& $OUTPUTF
    endif
    if (-e $STATSDIFFF) chmod g+rw $STATSDIFFF
    chgrp fsdev $SUBJECTS_DIR/ref.*
    chgrp fsdev $SUBJECTS_DIR/tst.*
    chmod g+rw $SUBJECTS_DIR/ref.*
    chmod g+rw $SUBJECTS_DIR/tst.*
  endif # ($RunIt)
end

set cmd=(asegstatsdiff ref_subj $TEST_SUBJ $LOG_DIR)
echo $cmd
set outfile=($LOG_DIR/asegstatsdiff.txt)
$cmd > $outfile
set stats_diff_status=$status
if ($stats_diff_status != 0) then
  printf "***FAILED :: asegstatsdiff (exit status=$stats_diff_status)\n" \
       >>& $OUTPUTF
  setenv FOUND_ERROR 1
  # continue running tests
else
  printf "   pass :: asegstatsdiff\n" >>& $OUTPUTF
endif

foreach hemi (rh lh)
    foreach parc (aparc aparc.a2009s)
        foreach meas (area volume thickness)
        set cmd=(aparcstatsdiff \
            ref_subj $TEST_SUBJ $hemi $parc $meas $LOG_DIR)
        echo $cmd
        set outfile=($LOG_DIR/aparcstatsdiff-$hemi-$parc-$meas.txt)
        $cmd > $outfile
        set stats_diff_status=$status
        if ($stats_diff_status != 0) then
            printf "***FAILED :: aparcstatsdiff-$hemi-$parc-$meas (exit status=$stats_diff_status)\n" \
                >>& $OUTPUTF
            setenv FOUND_ERROR 1
            # continue running tests
        else
            printf "   pass :: aparcstatsdiff-$hemi-$parc-$meas\n" >>& $OUTPUTF
        endif
        end
    end
end


#
# Organize the pile of results
#
cd $LOG_DIR
mkdir -p mri_diff > /dev/null
rm -Rf mri_diff/* > /dev/null
mkdir -p mris_diff > /dev/null
rm -Rf mris_diff/* > /dev/null
mkdir -p aseg > /dev/null
rm -Rf aseg/* > /dev/null
mkdir -p aparc > /dev/null
rm -Rf aparc/* > /dev/null

mv mri_diff* mri_diff/ > /dev/null
mv mris_diff* mris_diff/ > /dev/null
mv *seg* aseg/ > /dev/null
mv aparc* aparc/ > /dev/null
mv *.aparc* aparc/ > /dev/null

#
# move completed run to another dir, indicating that recon-all test
# has completed (beginning of this script checks if subject dir exists,
# and doesnt start recon-all if it does)
#
if ( ! $?SKIP_RECON) then
    rm -Rf $SUBJECTS_DIR/$TEST_SUBJ-done
    mv $SUBJECTS_DIR/$TEST_SUBJ $SUBJECTS_DIR/$TEST_SUBJ-done
    chgrp -R fsdev $SUBJECTS_DIR/$TEST_SUBJ-done
    chmod -R g+rw $SUBJECTS_DIR/$TEST_SUBJ-done
endif


#
# completion message (success or failure)
#

done:
chgrp -R fsdev $LOG_DIR
chmod -R g+rw $LOG_DIR
chgrp fsdev $OUTPUTF
chmod g+rw $OUTPUTF
set END_TIME=`date`
echo "Finish: $END_TIME" >>& $OUTPUTF
if ($?FOUND_ERROR) then
  echo "FAILURE(S) found in test_recon-all $TEST_SUBJ on $PROC (${HOST})" \
    >>& $OUTPUTF
  echo "Consult diff logs in $LOG_DIR for details" >>& $OUTPUTF
  echo "Output SUBJECTS_DIR is $SUBJECTS_DIR" >>& $OUTPUTF
  if ($RunIt) then
    mail -v -s "test_recon-all $TEST_SUBJ FAILURE(s) on $PROC (${HOST})" \
        $FMAIL_LIST < $OUTPUTF
    touch $SUBJECTS_DIR/test_recon-all_FAILED
    chmod a+w $SUBJECTS_DIR/test_recon-all_FAILED
    exit 1
  endif
else
  echo "Success running test_recon-all $TEST_SUBJ on $PROC (${HOST})" \
    >>& $OUTPUTF
  if ($RunIt) then
    mail -v -s "test_recon-all $TEST_SUBJ success on $PROC (${HOST})" \
      $MAIL_LIST < $OUTPUTF
    exit 0
  endif
endif
