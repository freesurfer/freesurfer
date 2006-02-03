#!/bin/tcsh -f

set vol_in = $1
set vol_out = $2

mri_convert ${vol_in} input.mnc
nu_correct -clobber input.mnc nu_1.mnc 
nu_correct -clobber nu_1.mnc nu_2.mnc 

mri_convert nu_2.mnc ${vol_out} 

rm -f input*
rm -f nu_1.mnc
rm -f *.imp
rm -f nu_2.*
