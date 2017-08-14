#!/bin/tcsh -ef

if (! $?base) then
  echo must set base to point to $SUBJECTS_DIR/$subject
  exit 1
endif

set odir = $base/mri/orig
set ldir = $base/label
set labels = (wm gm fluid)
set optdir = $odir/opt_files
mkdir -p $optdir

echo "set flist = $optdir/flist.dat"
set flist = $optdir/flist.dat
if (-e $flist) then
  rm $flist
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
    if ($fonly == opt.unmasked) then
	continue ;
    endif
    if ($fonly == fluid_mask) then
	continue ;
    endif
    if ($fonly == brainmask) then
	continue ;
    endif
    if ($fonly == sse) then
	continue ;
    endif
    if ($fonly == faf) then
	continue ;
    endif
    if ($fonly == bias) then
	continue ;
    endif
    
    echo "${f:r}  $flist"
    echo ${f:r} >> $flist
    foreach l ($labels)
            mri_label_vals $f $ldir/${l}.label > $odir/opt_files/$fonly.$l.dat
    end
end

