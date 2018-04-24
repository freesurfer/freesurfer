#!/bin/tcsh -ef
set s=FCD-JO
set hemi=rh
./mris_extract_patches -hemi $hemi -sdir $FCD
