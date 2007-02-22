#!/bin/tcsh -ef

set indir = /space/structure/1/freediffusion/data/practice/mri/dti/004
set infile = $indir/f.bhdr 
set outdir = /Users/ayendiki/fd-tmp/c
set outfmt = bhdr

set mask = /Users/ayendiki/fd-tmp/mask0_000.bfloat

set infofile = $indir/f-infodump.dat

set bvalue = `grep "sDiffusion\.alBValue\[1\]" $infofile | awk '{print $3}'`
set nt2    = `grep "sWiPMemBlock\.alFree\[8\]" $infofile | awk '{print $3}'`
set ndir = `grep "sDiffusion\.lDiffDirections" $infofile | awk '{print $3}'`
set diffmode = `grep "sWiPMemBlock\.alFree\[1\]" $infofile | awk '{print $3}'` 
switch($diffmode)
  case 3:
    set gradfile = /usr/local/freesurfer/dev/freediffusion/diffusion-toolbox/directions/6-cube-mghnew.txt
  breaksw
endsw

set cmd = (./dmri_tensoreig \
              --i $infile  \
              --m $mask \
              --b $bvalue --nacq $nt2 --ndir $ndir --g $gradfile \
              --o $outdir --ofmt $outfmt)
echo $cmd
$cmd

