#!/bin/tcsh -ef

if (! $?base) then
  echo must set base to point to $SUBJECTS_DIR/$subject
  exit 1
endif

set odir = $base/mri/orig
set ldir = $base/label
set labels = (wm gm fluid)
mkdir -p $odir/opt_files

if (-e flist.dat) then
  rm flist.dat
endif

set echo=1
foreach f ($odir/norm/*.mgz)
    set fonly = ${f:t}
    set fonly = ${fonly:r}
    if ($fonly == 001) then
	continue ;
    endif
    if ($fonly == opt) then
	continue ;
    endif
    if ($fonly == fluid_mask) then
	continue ;
    endif
    if ($fonly == brainmask) then
	continue ;
    endif
    if ($fonly == bias) then
	continue ;
    endif
    echo ${f:r} >> flist.dat
    foreach l ($labels)
            mri_label_vals $f $ldir/${l}.label > $odir/opt_files/$fonly.$l.dat
    end
end

