#!/bin/tcsh -f

#
# rebuild_gca_atlas.csh
#
# This script builds the subcortical atlas from a set of manually labelled
# training data sets.
# The user should call the script under a cluster machine (seychelles).
# In addition, a "subjects.csh" script should be available under the
# $SUBJECTS_DIR/scripts directory, which defines the training subject names
# and the first subject ($ONE_SUBJECT) to be used for building the initial
# template. Eg. /space/dijon/32/users/xhan/RBnew/scripts/subjects.csh
# A talairach registration should be generated for this first subject as well,
# in order to align the final atlas to the Talairach space.
# The registration file should be put under
# $SUBJECTS_DIR/$ONE_SUBJECT/mri/transforms/$TAL_MAN
# The script assumes the existence of following files under
# each subject's mri directory: nu.mgz, brain.mgz, nu_noneck.mgz and
# aseg_edited.mgz, where aseg_edited.mgz is the manually labelled volume
# from which the training data is derived.
# The final atlas is stored as $SUBJECTS_DIR/average/RB_all_`date +%F`.gca
# "nu_noneck.mgz" is needed for every subject to build the gca with skull.
#
# Original author: Xiao Han
# CVS Revision Info:
#    $Author: nicks $
#    $Date: 2007/04/25 18:52:23 $
#    $Revision: 1.10 $
#
# Copyright (C) 2002-2007,
# The General Hospital Corporation (Boston, MA).
# All rights reserved.
#
# Distribution, usage and copying of this software is covered under the
# terms found in the License Agreement file named 'COPYING' found in the
# FreeSurfer source code root directory, and duplicated here:
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
#
# General inquiries: freesurfer@nmr.mgh.harvard.edu
# Bug reports: analysis-bugs@nmr.mgh.harvard.edu
#


set VERSION='$Id: rebuild_gca_atlas.csh,v 1.10 2007/04/25 18:52:23 nicks Exp $';

#set echo=1

# these are the subjects to use in training.
# one of them, the target, set to ONE_SUBJECT, should be declared is this file:
source ${SUBJECTS_DIR}/scripts/subjects.csh

# created GCA files, and log file, get this date appended
set DATE=(`date +%F`)

# pbs special configure options:
set PBCONF="-l nodes=1:opteron -m $USER"

# optionally choose to not run commands (but commands are echoed
# in logfile) by setting RunIt=0
set RunIt=1

# optionally jump to creation of input volumes (skipping atlas build)
if ("$1" == "create_inputs") goto create_inputs

##########################################################################
#
# Log file:
#
set LF=(${SUBJECTS_DIR}/rebuild_gca_atlas_${DATE}.log)
if ("$1" != "2") then
    echo "Log file is $LF"
    echo "\n\n========= rebuild_gca_atlas ==============\n" >& $LF
    echo "\n$VERSION\n" >>& $LF
    echo "Start: `date`\n" >>& $LF
    if ( ! $RunIt) echo "Command execution is OFF (RunIt=0)"
endif
##########################################################################



##########################################################################
#
# Binaries needed:
#
set emreg=(mri_em_register)
set careg=(mri_ca_register)
set canorm=(mri_ca_normalize)
set train=(mri_ca_train)
$emreg  --all-info >>& $LF
$careg  --all-info >>& $LF
$canorm --all-info >>& $LF
$train  --all-info >>& $LF
echo "\n\n" >>& $LF
##########################################################################



##########################################################################
#
# Inputs:
#
set SEG_VOL=(aseg_edited.mgz) # filename for manual segmentation
set ORIG_VOL=(nu.mgz)
set MASK_VOL=(brain.mgz) # filename for brain mask
set T1_NONECK=(nu_noneck.mgz) # file to build the atlas gca_with_skull
set TAL_MAN=(RB_AVERAGE3new.xfm) # optional manual tal registration file
set INPUTS=(${SEG_VOL} ${ORIG_VOL} ${MASK_VOL} ${T1_NONECK})
set ALL_SUBJS=(${SUBJECTS} ${ONE_SUBJECT})
foreach subject (${ALL_SUBJS}) # check for existence of required inputs
    foreach input (${INPUTS})
        if ( ! -e ${SUBJECTS_DIR}/$subject/mri/$input ) then
            echo "Missing ${SUBJECTS_DIR}/$subject/mri/$input!"
            exit 1
        endif
    end
end
if ( ! -e ${SUBJECTS_DIR}/${ONE_SUBJECT}/mri/transforms/${TAL_MAN} ) then
    echo "Missing ${SUBJECTS_DIR}/${ONE_SUBJECT}/mri/transforms/${TAL_MAN}!"
    exit 1
endif

##########################################################################



##########################################################################
#
# Outputs:
#
mkdir -p ${SUBJECTS_DIR}/average
set GCA=(${SUBJECTS_DIR}/average/RB_all_${DATE}.gca)
set GCA_ONE=(${SUBJECTS_DIR}/average/RB_one_${DATE}.gca)
set GCA_SKULL=(${SUBJECTS_DIR}/average/RB_all_withskull_${DATE}.gca)

set LTA_ONE=(talairach_one.lta)
set M3D_ONE=(talairach_one.m3z)

set LTA=(talairach.lta)
set M3D=(talairach.m3z)

set T1_VOL = norm.mgz
##########################################################################



# optionally divided into two stages.
# if the argument to this script is '2', then it skips the first stage:
if ("$1" == "2") goto stage2



##########################################################################
#
# Normalize brains:
#
echo "$canorm each subject using its ${SEG_VOL}, producing ${T1_VOL}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    if(-e $mridir/${T1_VOL}) rm -f $mridir/${T1_VOL}

    set cmd=($canorm -mask $mridir/${MASK_VOL})
    set cmd=($cmd -seg $mridir/${SEG_VOL})
    set cmd=($cmd $mridir/${ORIG_VOL} noatlas noxform $mridir/${T1_VOL})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until $canorm has finished each subject
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/${T1_VOL}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $canorm on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# Train using ONE_SUBJECT
#
if( -e ${GCA_ONE} ) rm -f ${GCA_ONE}
if(-e ${SUBJECTS_DIR}/$ONE_SUBJECT/mri/transforms/${TAL_MAN}) then
    echo "Using manual talairach registration ${TAL_MAN}..." >>& $LF
    set MAN_TAL=(-xform ${TAL_MAN})
else
    echo "Initial talairach registration unavailable (no xform)..." >>& $LF
    set MAN_TAL=
endif
echo "$train using one subject: $ONE_SUBJECT, producing ${GCA_ONE}..."
set cmd=($train -prior_spacing 2 -node_spacing 8 -mask ${MASK_VOL})
set cmd=($cmd -parc_dir ${SEG_VOL} ${MAN_TAL} -T1 ${T1_VOL})
set cmd=($cmd $ONE_SUBJECT ${GCA_ONE})
echo $cmd >>& $LF
if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
echo "\n\n" >>& $LF
# waiting
set TEST = 0
if ($RunIt) set TEST = 1
while($TEST)
    if(-e ${GCA_ONE}) set TEST = 0
    sleep 30
end
echo "\t...finished $train, produced $GCA_ONE" >>& $LF
##########################################################################



##########################################################################
#
# EM_registration by GCA_ONE
#
echo "$emreg each subject to one-subj GCA, producing transforms/${LTA_ONE}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    mkdir -p $mridir/transforms
    set tdir=($mridir/transforms)
    if ( -e $tdir/${LTA_ONE} ) rm -f $tdir/${LTA_ONE}

    set cmd=($emreg -mask $mridir/${MASK_VOL} $mridir/${ORIG_VOL})
    set cmd=($cmd ${GCA_ONE} $mridir/transforms/${LTA_ONE})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the registration is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/transforms/${LTA_ONE}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $emreg on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# normalization by GCA_ONE
#
echo "$canorm each subject using one-subject GCA, producing ${T1_VOL}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    if ( -e $mridir/${T1_VOL} ) rm -f $mridir/${T1_VOL}

    set cmd=($canorm -mask $mridir/${MASK_VOL} $mridir/${ORIG_VOL})
    set cmd=($cmd ${GCA_ONE} $mridir/transforms/${LTA_ONE})
    set cmd=($cmd $mridir/${T1_VOL})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the normalization is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/${T1_VOL}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $canorm on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# CA_registration by GCA_ONE
#
echo "$careg each subject to one-subj GCA, producing transforms/${M3D_ONE}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    if ( -e  $mridir/transforms/${M3D_ONE}) rm -f $mridir/transforms/${M3D_ONE}

    set cmd=($careg -smooth 1.0 -levels 2 -mask $mridir/${MASK_VOL})
    set cmd=($cmd -T $mridir/transforms/${LTA_ONE} $mridir/${T1_VOL})
    set cmd=($cmd ${GCA_ONE} $mridir/transforms/${M3D_ONE})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the registration is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/transforms/${M3D_ONE}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $careg on subject $subject" >>& $LF
end
##########################################################################


#
# optionally stop here, and adjust SUBJECTS to contain
# only well-aligned subjects.
# to stop here, use '1' as the input argument to this script.
# to start-up here, use '2' as the input arg. default is run both stages.
#
if ("$1" == "1") exit 0
stage2:


##########################################################################
#
# TRAIN FROM SEGMENTED_SUBJECTS USING M3D_ONE
#
if (-e ${GCA} ) rm -f ${GCA}

echo "$train using all subjects, using ${M3D_ONE}, producing ${GCA}..."
set cmd=($train -prior_spacing 2 -node_spacing 4 -mask ${MASK_VOL})
set cmd=($cmd -parc_dir ${SEG_VOL} -xform ${M3D_ONE} -T1 ${T1_VOL})
set cmd=($cmd ${SUBJECTS} ${GCA})
echo $cmd >>& $LF
if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
echo "\n\n" >>& $LF
# waiting
set TEST = 0
if ($RunIt) set TEST = 1
while($TEST)
    if(-e ${GCA}) set TEST = 0
    sleep 30
end
echo "\t...finished $train, produced ${GCA}" >>& $LF
##########################################################################



##########################################################################
#
# REGISTER ALL BRAINS TO GCA
# EM_registration by GCA
#
echo "$emreg each subject to GCA, producing transforms/${LTA}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    set tdir=($mridir/transforms)
    if ( -e $tdir/${LTA} ) rm -f $tdir/${LTA}

    set cmd=($emreg -mask $mridir/${MASK_VOL} $mridir/${ORIG_VOL})
    set cmd=($cmd ${GCA})
    set cmd=($cmd $mridir/transforms/${LTA})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the registration is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/transforms/${LTA}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $emreg on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# normalization by GCA
#
echo "$canorm each subject using GCA, producing ${T1_VOL}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    if ( -e $mridir/${T1_VOL} ) rm -f $mridir/${T1_VOL}

    set cmd=($canorm -mask $mridir/${MASK_VOL} $mridir/${ORIG_VOL})
    set cmd=($cmd ${GCA} $mridir/transforms/${LTA})
    set cmd=($cmd $mridir/${T1_VOL})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the normalization is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/${T1_VOL}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $canorm on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# CA_registration by GCA
#
echo "$careg each subject to GCA, producing transforms/${M3D}..."
foreach subject (${SUBJECTS})
    set mridir=(${SUBJECTS_DIR}/$subject/mri)
    if ( -e $mridir/transforms/${M3D} ) rm -f $mridir/transforms/${M3D}

    set cmd=($careg -smooth 1.0 -mask $mridir/${MASK_VOL})
    set cmd=($cmd -T $mridir/transforms/${LTA} $mridir/${T1_VOL})
    set cmd=($cmd ${GCA} $mridir/transforms/${M3D})
    echo $cmd >>& $LF
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
echo "\n\n" >>& $LF
# wait until the registration is done
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/transforms/${M3D}) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished $careg on subject $subject" >>& $LF
end
##########################################################################



##########################################################################
#
# RETRAIN GCA USING M3D
#
rm -f ${GCA}

echo "$train, using ${M3D}, producing ${GCA}..."
set cmd=($train -prior_spacing 2 -node_spacing 4 -mask ${MASK_VOL})
set cmd=($cmd -parc_dir ${SEG_VOL} -xform ${M3D} -T1 ${T1_VOL})
set cmd=($cmd ${SUBJECTS} ${GCA})
echo $cmd >>& $LF
if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
echo "\n\n" >>& $LF
# waiting
set TEST = 0
if ($RunIt) set TEST = 1
while($TEST)
    if(-e ${GCA}) set TEST = 0
    sleep 30
end
echo "\t...finished $train, produced final $GCA" >>& $LF

rm -f ${GCA_SKULL}

echo "$train, using ${LTA}, producing ${GCA_SKULL}..."
set cmd=($train -prior_spacing 2 -node_spacing 4)
set cmd=($cmd -parc_dir ${SEG_VOL} -xform ${LTA} -T1 ${T1_NONECK})
set cmd=($cmd ${SUBJECTS} ${GCA_SKULL})
echo $cmd >>& $LF
if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
echo "\n\n" >>& $LF
# waiting
set TEST = 0
if ($RunIt) set TEST = 1
while($TEST)
    if(-e ${GCA_SKULL}) set TEST = 0
    sleep 30
end
echo "\t...finished $train, produced final $GCA_SKULL" >>& $LF
##########################################################################

echo "Finished `date`" >>& $LF

exit 0




##########################################################################
#
# Run this if inputs do not exist or need updating.
#
# Usage:  rebuild_gca_atlas.csh create_inputs
#

create_inputs:

echo "create <orig.mgz>, nu.mgz, nu_noneck.mgz and brain.mgz..."
foreach subject (${SUBJECTS}) # check for existence of orig.mgz
    if ( ! -e ${SUBJECTS_DIR}/$subject/mri/orig.mgz ) then
        if ( -e ${SUBJECTS_DIR}/$subject/mri/orig.img ) then
            set cmd=(mri_convert ${SUBJECTS_DIR}/$subject/mri/orig.img)
            set cmd=($cmd ${SUBJECTS_DIR}/$subject/mri/orig.mgz)
            if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
        else
            echo "Missing ${SUBJECTS_DIR}/$subject/mri/orig*"
            exit 1
        endif
    endif
end
# wait until finished each subject
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/orig.mgz) then
            set TEST = 0
        else
            sleep 1
        endif
    end
    echo "\t...subject $subject has orig.mgz"
end
foreach subject (${SUBJECTS})
    rm -f ${SUBJECTS_DIR}/$subject/mri/brain.mgz
    set cmd=(recon-all -s $subject)
    set cmd=($cmd -nuintensitycor)
    set cmd=($cmd -talairach)
    set cmd=($cmd -normalization)
    set cmd=($cmd -skullstrip)
    set cmd=($cmd -subcortseg)
    set cmd=($cmd -normalization2)
    if ($RunIt) pbsubmit ${PBCONF} -c "$cmd"
end
# wait until finished each subject
foreach subject (${SUBJECTS})
    set TEST = 0
    if ($RunIt) set TEST = 1
    while($TEST)
        if(-e ${SUBJECTS_DIR}/$subject/mri/brain.mgz) then
            set TEST = 0
        else
            sleep 30
        endif
    end
    echo "\t...finished subject $subject"
end
##########################################################################

