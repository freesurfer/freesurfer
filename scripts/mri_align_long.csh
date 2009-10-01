#! /bin/tcsh -f

#
# mri_align_long.csh
#
# aligns all tps norm.mgz and aseg.mgz to the base space
#
# Original Author: Martin Reuter
# CVS Revision Info:
#    $Author: mreuter $
#    $Date: 2009/10/01 15:27:30 $
#    $Revision: 1.1 $
#
# Copyright (C) 2005-2008,
# The General Hospital Corporation (Boston, MA).
# All rights reserved.
#
# Distribution, usage and copying of this software is covered under the
# terms found in the License Agreement file named 'COPYING' found in the
# FreeSurfer source code root directory, and duplicated here:
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
#
# General inquiries and bug reports: freesurfer@nmr.mgh.harvard.edu
#


set VERSION = '$Id: mri_align_long.csh,v 1.1 2009/10/01 15:27:30 mreuter Exp $';


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

set baseid  = $1
set basedir = $SUBJECTS_DIR/$baseid
set tps = `cat $basedir/base-tps`

foreach s ($tps)

  echo TP: $s.long.$baseid

  set longdir = $SUBJECTS_DIR/$s.long.$baseid
	set lta = $basedir/mri/transforms/${s}_to_${baseid}.lta

	set src = $longdir/mri/norm.mgz
	set trg = $longdir/mri/norm-base.mgz
	set cmd = (mri_convert -at $lta -rl $basedir/mri/norm.mgz $src $trg)
  echo "\n $cmd \n"
  eval $cmd

	set src = $longdir/mri/aseg.mgz
	set trg = $longdir/mri/aseg-base.mgz
	set cmd = (mri_convert -at $lta -rl $basedir/mri/norm.mgz -rt nearest \
	                $src $trg)
  echo "\n $cmd \n"
  eval $cmd

end
