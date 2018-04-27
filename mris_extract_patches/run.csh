#!/bin/tcsh -ef
source $FCD2/scripts/subjects.csh
set echo=1

set sd=$SUBJECTS_DIR
foreach sbase ($SURF_AUTO_LABEL)
    set s=${sbase}.gentle
    if (-e $SUBJECTS_DIR/$s/label/rh.${FCD_LABEL}.label) then
	set hemi=rh
    else if (-e $SUBJECTS_DIR/$s/label/lh.${FCD_LABEL}.label) then
	set hemi=lh
    else
	exit 1
    endif
    
    set odir=$sd/$s/patches
    if (-e $odir == 0) then
	mkdir $odir
    endif
    
    ./mris_extract_patches -w 16 -l ${FCD_LABEL} -hemi $hemi -sdir $sd $s  $odir&

    echo vglrun freeview -v $sd/$s/mri/norm.mgz  -f $sd/$s/surf/$hemi.white:label=$sd/$s/label/$hemi.${FCD_LABEL}.label $sd/$s/surf/$hemi.inflated:label=$sd/$s/label/$hemi.${FCD_LABEL}.label&
    echo vglrun freeview -v $odir/$hemi.patches.mgz&
    sleep 30
#    break
end


