#!/bin/tcsh -ef
# usage: make_exvivo_surfs.csh <subject name> <input samseg> <input intensity vol>

if ( $#argv < 4) then
    usage: make_exvivo_surfs.csh <subject name> <input samseg> <input intensity vol>
    exit 1
endif

set s=$1
set samseg=$2
set norm=$3
set bdir=$SUBJECTS_DIR/$s
set mdir=$bdir/mri
set sdir=$bdir/surf

if (-e $bdir == 0) then
    echo creating subject directory tree
    mksubjdirs $s
endif



#mri_mask -dilate 3 -no_cerebellum -lh $mdir/flash20_polyfit2_masked_v2.RegOCT.mgz $mdir/flash20_polyfit2_masked_v2_600.RegOCT_masked_crispSegmentation.nii $mdir/flash20.150um.RegOCT.masked.dilated.mgz
mri_convert $samseg $mdir/aseg.mgz
mri_mask -no_cerebellum -lh $norm $samseg $mdir/norm.mgz
mri_extract_largest_CC -l 2 $samseg $mdir/wm.seg.mgz
mri_binarize --binval 255 --i $mdir/wm.seg.mgz --o $mdir/wm.seg.bin.mgz --min 1 --max 3
mri_convert -odt uchar $mdir/wm.seg.bin.mgz $mdir/wm.seg.bin.mgz
mri_edit_wm_with_aseg -lh -keep-in $mdir/wm.seg.bin.mgz $norm $mdir/aseg.mgz $mdir/wm.asegedit.mgz 
mri_extract_largest_CC -l 255 $mdir/wm.asegedit.mgz $mdir/wm.asegedit.cc.mgz
mri_pretess $mdir/wm.asegedit.mgz wm $mdir/norm.mgz $mdir/wm.mgz 
mri_fill -lhonly -a ../scripts/ponscc.cut.log -segmentation $mdir/aseg.mgz $mdir/wm.mgz $mdir/filled.mgz 
mri_tessellate $mdir/filled.mgz 255 $sdir/lh.orig
#mri_mc $mdir/filled.mgz 255 $sdir/lh.orig
mris_smooth $sdir/lh.orig $sdir/lh.smoothwm
