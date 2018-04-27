#!/bin/tcsh -ef
source $FCD2/scripts/subjects.csh
set echo=1

set sd=$SUBJECTS_DIR
set max_running=15
foreach sbase ($SURF_AUTO_LABEL)
    set s=${sbase}.gentle
    set hemi=lh
    
    set odir=$sd/$s/cortex_patches
    if (-e $odir == 0) then
	mkdir $odir
    endif

    set bin=./mris_extract_patches
    set nrunning = (`ps aux | grep $bin | wc -l`)
    while ($nrunning > $max_running)
	echo $nrunning running - sleeping
	sleep 60
	set nrunning = (`ps aux | grep $bin | wc -l`)
    end
    $bin -w 16 -l cortex -hemi $hemi -sdir $sd $s  $odir&

    echo vglrun freeview -v $sd/$s/mri/norm.mgz  -f $sd/$s/surf/$hemi.white:label=$sd/$s/label/$hemi.cortex.label $sd/$s/surf/$hemi.inflated:label=$sd/$s/label/$hemi.cortex.label
    echo vglrun freeview -v $odir/$hemi.patches.mgz&
end


