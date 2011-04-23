#!/bin/bash
#
# Set up all the links to the data that we need
#
# If this script failed or was stopped, do 
# 
#  rm -r */
# 
# to clean up the mess and start with a clean plate
#


# Define a function that will show what we're doing, and
# exits if something goes wrong
function doIt {
  echo $1
  eval $1
  if [ $? != 0 ]; then
    echo "failed to do $1"
  exit -1
fi

}


#
# The following hard-coded stuff was compiled from the descriptions available at http://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalLabeling,
# in addition to the knowledge that all data is sym-linked to /autofs/space/rock_004/users/hires, and that the subdirectory "recons" there
# contains in vivo scans (directory names _invivo, _InVivo, etc) (thanks Allison!), as well as low-res ex vivo scans (information obtained
# from http://surfer.nmr.mgh.harvard.edu/fswiki/ExvivoRecons)
#
# We name our subjects in the same order as in the "Labeling Log" table. Furthermore, if more than one manual segmentation is available, we
# number the raters in the order in which they appear in the same table as well.


# -----------------------------------------------------
# FHS18 -- subject01
# -----------------------------------------------------

subjectName="subject01"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS18/FHS18_ec_031408/mri/flash/flash20_EC_100um_avg.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS18/FHS18_ec_031408/mri/flash/HP_06152009LC.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS18_Invivo/mri"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
# I9 -- subject02
# -----------------------------------------------------

subjectName="subject02"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I9/I9_mtl_may_3_2004/mri/parameter_maps/I9/flash23.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I9/I9_mtl_may_3_2004/mri/parameter_maps/I9/HP_07072009LC.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo for this subject, so nothing to do here

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
# FHS20 -- subject03
# -----------------------------------------------------

subjectName="subject03"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS20/FHS20_lh_ec_120um/mri/flash/flash20_ec_120um_avg.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS20/FHS20_lh_ec_120um/mri/flash/HP_06102009NR.mgz seg_rater1.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS20/FHS20_lh_ec_120um/mri/flash/HP_07012009MK.mgz seg_rater2.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS20_invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
# FHS8 -- subject04
# -----------------------------------------------------

subjectName="subject04"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS8_RH/FHS8_rh_150um_04_09_2007/mri/flash/flash_avg_ras.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS8_RH/FHS8_rh_150um_04_09_2007/mri/flash/HP_06252009NR.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS8_Invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
# FHS21 -- subject05
# -----------------------------------------------------

subjectName="subject05"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS21/FHS21_hippo_amyg/mri/flash/flash20_ec_amyg_120um_avg.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS21/FHS21_hippo_amyg/mri/flash/HP_06102009RK.mgz seg_rater1.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS21/FHS21_hippo_amyg/mri/flash/HP_06262009LC.mgz seg_rater2.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS21_Invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
# FHS5 -- subject06
# -----------------------------------------------------

subjectName="subject06"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS5/FHS5_rh_ec_150um_11_6_2006/mri/flash/flash20_150um_EC_avg_ras.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS5/FHS5_rh_ec_150um_11_6_2006/mri/flash/HP_06242009RK.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid 
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS5_in_vivo_recon/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
# FHS19 -- subject07
# -----------------------------------------------------

subjectName="subject07"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS19/FHS19_ec_amyg_032808/mri/flash/flash20_ec_amyg_100um_avg.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS19/FHS19_ec_amyg_032808/mri/flash/HP_06232009CP.mgz seg_rater1.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS19/FHS19_ec_amyg_032808/mri/flash/HP_0818_sk.mgz seg_rater2.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid 
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS19_invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
# I13 rh -- subject08
# -----------------------------------------------------

subjectName="subject08"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I13/I13_whole_hemi_hippo_4_1_2005/mri/flash20_250um_avg_rotated_jca_07292009.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I13/I13_whole_hemi_hippo_4_1_2005/mri/HP_07302009CP_rh.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/I13_BWH_in_vivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
# I13 lh -- subject09
# -----------------------------------------------------

subjectName="subject09"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I13/I13_whole_hemi_hippo_4_1_2005/mri/flash20_250um_avg_rotated_cmp_08102009.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I13/I13_whole_hemi_hippo_4_1_2005/mri/HP_08102009CP_lh.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid 
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/I13_BWH_in_vivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."



# -----------------------------------------------------
#  I21 -- subject10
# -----------------------------------------------------

subjectName="subject10"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I21/I21_lh_mtl_150um_me2_09_22_2006/mri/flash/mef20_150um_avg_rotated_cmp08132009.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I21/I21_lh_mtl_150um_me2_09_22_2006/mri/flash/HP_08132009CP.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/I21_lh_recon/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."


# -----------------------------------------------------
#  I28 -- subject11
# -----------------------------------------------------

subjectName="subject11"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I28/I28_rh_ec_amyg_101708/mri/flash/flash20_ec_100um_avg.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I28/I28_rh_ec_amyg_101708/mri/flash/HP_07282009MK.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."


# -----------------------------------------------------
#  I23 -- subject12
# -----------------------------------------------------

subjectName="subject12"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I23/I23_lh_mtl_150um_11_27_2006/mri/flash/flash20_MTL_150um_avg_oriented_081809KH.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I23/I23_lh_mtl_150um_11_27_2006/mri/flash/HP_12292009NR.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I10 -- subject13: TODO this one is not correctly symlinked to /autofs/space/rock_004/users/hires
# -----------------------------------------------------

subjectName="subject13"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/cluster/con_006/users/hires/I10_9_21_2004_whole_hemi_150um_labeled_volume/flash_low_BW_avg_oriented_SK.mgz orig.mgz"
doIt "ln -s /autofs/cluster/con_006/users/hires/I10_9_21_2004_whole_hemi_150um_labeled_volume/HP_12292009MR.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I27 -- subject14
# -----------------------------------------------------

subjectName="subject14"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I27/I27_mtl_small_sol/mri/flash/flash20_100um_avg_rotated_KH102809.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I27/I27_mtl_small_sol/mri/flash/HP_01192010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  FHS34 -- subject15: TODO this one is not correctly symlinked to /autofs/space/rock_004/users/hires
# -----------------------------------------------------

subjectName="subject15"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/cluster/exvivo/FHS34_lh_hp_150um_temp/flash20_hp_150um_avg_ras_rotated_hp_04152010MR.mgz orig.mgz"
doIt "ln -s /autofs/cluster/exvivo/FHS34_lh_hp_150um_temp/HP_06282010NR.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS34_Invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I26 -- subject16
# -----------------------------------------------------

subjectName="subject16"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I26//I26_rh_EC_07_06_2007/mri/flash/HP_04262010_KN/flash20_120um_ds_150_ras_rot_cropped.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I26//I26_rh_EC_07_06_2007/mri/flash/HP_04262010_KN/HP_04262010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I_ellen -- subject17
# -----------------------------------------------------

subjectName="subject17"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I_ellen/I_ellen_aug_16_04_mtl/mri/flash200um_avg_ras_oriented_6152010KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I_ellen/I_ellen_aug_16_04_mtl/mri/HP_06152010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  FHS38 -- subject18
# -----------------------------------------------------

subjectName="subject18"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS38/FHS38_lh_mtl_120um/mri/flash/flash20_mtl_120um_avg_ds_150um_ras_rot_crop_AS.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS38/FHS38_lh_mtl_120um/mri/flash/HP_06282010CP.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS38_Invivo/mri/"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I25 -- subject19
# -----------------------------------------------------

subjectName="subject19"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I25/I25_lh_ec_06_01_2007/mri/flash/flash20_EC_120um_avg_ras_rotated_hp_6302010MR.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I25/I25_lh_ec_06_01_2007/mri/flash/HP_06302010MR.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I16 -- subject20
# -----------------------------------------------------

subjectName="subject20"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I16/I16_lh_150um_EC_08_24_06/mri/flash/flash20_MTL_150um_avg_ras_rot_crop_KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I16/I16_lh_150um_EC_08_24_06/mri/flash/HP_7092010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  I18 -- subject21
# -----------------------------------------------------

subjectName="subject21"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I18/I18_jun_6_whole_hemi_mtl_150um/mri/flash20_150um_avg_ras_rot_crop_KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I18/I18_jun_6_whole_hemi_mtl_150um/mri/HP_07122010NR.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  FHS37 -- subject22
# -----------------------------------------------------

subjectName="subject22"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS37/FHS37_rh_mtl_120um_rescan/mri/flash/flash20_mtl_120um_avg_ras_rot_crop2_KN_CP.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS37/FHS37_rh_mtl_120um_rescan/mri/flash/HP_07152010CP.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."




# -----------------------------------------------------
#  FHS39 -- subject23
# -----------------------------------------------------

subjectName="subject23"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS39/FHS39_lh_mtl_120um/mri/flash/flash20_mtl_120um_avg_ras_crop_KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS39/FHS39_lh_mtl_120um/mri/flash/HP_07212010NR.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir="/autofs/space/rock_004/users/hires/recons/FHS39_Invivo/mri"
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."





# -----------------------------------------------------
#  I29 -- subject24  
# -----------------------------------------------------

subjectName="subject24"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/I29/I29_RT_100um/mri/flash/flash20_100um_avg_ds_150um_rot_crop_KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/I29/I29_RT_100um/mri/flash/HP_8122010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."


# -----------------------------------------------------
#  FHS36_rh -- subject25  
# -----------------------------------------------------

subjectName="subject25"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS36/FHS36_rh_mtl_120um/mri/flash/flash20_mtl_120um_avg_ras_rot_crop_9082010KN.mgz orig.mgz"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS36/FHS36_rh_mtl_120um/mri/flash/HP_9082010KN.mgz seg_rater1.mgz"
doIt "cd .."

# No in vivo data for this subject

# Go back onto main directory
doIt "cd .."


# -----------------------------------------------------
#  FHS31_rh -- subject26: TODO the manual segmnentation HP_9172010MR.mgz does not exist  
# -----------------------------------------------------

subjectName="subject26"

# Make subject directory and move in it
doIt "mkdir $subjectName"
doIt "cd $subjectName"

# Make subdirectory for highres ex vivo and copy data and manual segmentations into it
doIt "mkdir exVivoHighRes"
doIt "cd exVivoHighRes"
doIt "ln -s /autofs/space/rock_004/users/hires/FHS31/FHS31_rh_mtl_120um/mri/flash/flash20_mtl_120um_avg_rotated_crop_hp_9172010MR.mgz orig.mgz"
doIt "ln -s /autofs/space/sake_021/users/hires/FHS31_rh_mtl_120um/mri/flash/HP_9172010MR_crop.mgz seg_rater1.mgz"
doIt "cd .."

# Make subdirectory for in vivo and copy the following data into it:
#   mri/rawavg.mgz  : the basic data without any preprocessing
#   mri/orig.mgz    : resampled to 1x1x1mm grid
#   mri/nu.mgz      : bias field corrected with N3
#   mri/transforms/talairach.lta : affine transform mapping the data onto FreeSurfer's mni305.cor.mgz template
#                                  (see Atlas3D/scripts/WholeBrainParcellation/README_CoregistrationWithAtlas.txt
#                                   for information on how to use that)
doIt "mkdir inVivo"
doIt "cd inVivo"
sourceDir=""
doIt "ln -s $sourceDir/rawavg.mgz ."
doIt "ln -s $sourceDir/orig.mgz ."
doIt "ln -s $sourceDir/nu.mgz ."
doIt "ln -s $sourceDir/transforms/talairach.lta ."
doIt "cd .."

# Go back onto main directory
doIt "cd .."

