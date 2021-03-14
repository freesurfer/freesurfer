#! /bin/tcsh -f

#
# mri_align_long.csh
#
# aligns all tps norm.mgz and aseg.mgz to the base space
#
# Original Author: Martin Reuter
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


set VERSION = 'mri_align_long.csh @FS_VERSION@';


if ($#argv < 1) then
   echo
   echo " mri_align_long.csh base-id "
	 echo
	 echo "   Aligns all longitudinal norm and aseg files to the base space."
	 echo "   Creates norm-base.mgz and aseg.base.mgz in each tpN.long.base/mri dir."
	 echo
	 echo
   exit 1;   
endif

echo Using SUBJECTS_DIR $SUBJECTS_DIR

set baseid  = $1
set basedir = $SUBJECTS_DIR/$baseid

if (! -e $basedir) then
  echo "$basedir does not exist!"
	echo "Please set correct SUBJECTS_DIR environment variable ..."
	exit 1
endif
if (! -e $basedir/base-tps) then
  echo "$basedir/base-tps does not exist!"
	echo "Is this really the base subject?"
	exit 1
endif
set tps = `cat $basedir/base-tps`


@ n = 0 ;
foreach s ($tps)
  @   n = $n + 1 ;

  echo
  echo TP: $s.long.$baseid

  set longdir = $SUBJECTS_DIR/$s.long.$baseid
	set lta = $basedir/mri/transforms/${s}_to_${baseid}.lta

	set src = $longdir/mri/norm.mgz
	if (! -e $src) then
	  echo "$src does not exist"
		exit 1
	endif
	set trg = $longdir/mri/norm-base.mgz
	set cmd = (mri_convert -at $lta -rl $basedir/mri/norm.mgz $src $trg)
  echo "\n $cmd \n"
  eval $cmd

	set src = $longdir/mri/aseg.mgz
	if (! -e $src) then
	  echo "$src does not exist"
		exit 1
	endif
	set trg = $longdir/mri/aseg-base.mgz
	set cmd = (mri_convert -at $lta -rl $basedir/mri/norm.mgz -rt nearest \
	                $src $trg)
  echo "\n $cmd \n"
  eval $cmd

end

echo The aligned norm-base.mgz and aseg-base.mgz can now be found in
echo "<tpN>.long.$baseid/mri"
echo

if ($n == 2) then
   echo Inspect with e.g.:
	 echo tkmedit -f $SUBJECTS_DIR/$tps[1].long.$baseid/mri/norm-base.mgz \
	              -aux $SUBJECTS_DIR/$tps[2].long.$baseid/mri/norm-base.mgz \
								-segmentation $SUBJECTS_DIR/$tps[1].long.$baseid/mri/aseg-base.mgz \
								-aux-segmentation $SUBJECTS_DIR/$tps[2].long.$baseid/mri/aseg-base.mgz
endif
