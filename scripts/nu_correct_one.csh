#!/bin/tcsh -ef


set s=$1
set iter=3

set mdir=$SUBJECTS_DIR/$s/mri
if (-e $mdir/orig == 0) then
		set ORIG = orig.mgz
else
		set ORIG = orig
endif

mri_convert $mdir/$ORIG /tmp/nu$$0.mnc
nu_correct -stop .0001 -iterations $iter -normalize_field -clobber /tmp/nu$$0.mnc /tmp/nu$$1.mnc
mkdir -p $mdir/nu
mri_convert /tmp/nu$$1.mnc $mdir/nu
rm /tmp/nu$$0.mnc /tmp/nu$$1.mnc
find  /tmp  -prune -name "*.imp" -user $LOGNAME -exec rm -f {} \;
#rm -f /tmp/*.imp


