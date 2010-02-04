#!/usr/bin/tcsh -ef
set echo

# CVS Info : $Id: test_fs_tutorial.csh,v 1.1 2010/02/04 19:37:09 krish Exp $
# This script executes all the command-line entries in all the tutorials of the FsTutorial
# It excludes all the tkmedit, tksurfer, qdec commands until the end where it invokes the 
# three commands in succession.
# Also excludes commands which take long ( more than a minute ) in the tutorial

# Notes to the user :
# Set FREESURFER_HOME before executing this script
# Also make sure $FREESURFER_HOME/subjects/buckner_data exists
# or is symlinked to the buckner data directory

#------------------------------
# ---- OutputData ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs

setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $SUBJECTS_DIR

#tkmedit good_output brainmask.mgz lh.white \
#        -aux T1.mgz -aux-surface rh.white \
#        -segmentation aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt

#tksurfer good_output lh inflated

#------------------------------
# ---- GroupAnalysis ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA/group_analysis_tutorial
cd $SUBJECTS_DIR/glm


mris_preproc --fsgd gender_age.fsgd \
  --cache-in thickness.fwhm10.fsaverage \
  --target fsaverage --hemi lh \
  --out lh.gender_age.thickness.10.mgh


mri_glmfit \
  --y lh.gender_age.thickness.10.mgh \
  --fsgd gender_age.fsgd dods\
  --C lh-Avg-thickness-age-Cor.mtx \
  --surf fsaverage lh \
  --cortex \
  --glmdir lh.gender_age.glmdir


#tksurfer fsaverage lh inflated \
#  -annot aparc.annot -fthresh 2 \
#-overlay lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/sig.mgh


mri_glmfit-sim \
  --glmdir lh.gender_age.glmdir \
  --sim mc-z 5 4 mc-z.negative \
  --sim-sign neg \
  --overwrite

cat lh.gender_age.glmdir/csd/mc-z.neg4.j001-lh-Avg-thickness-age-Cor.csd

cat lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/mc-z.neg4.sig.cluster.summary

#tksurfer fsaverage lh inflated \
#  -annot lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/mc-z.neg4.sig.ocn.annot \
#  -fthresh 2 -curv -gray

#------------------------------
# ---- QdecGroupAnalysis ----
#------------------------------


setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs

setenv SUBJECTS_DIR $TUTORIAL_DATA/group_analysis_tutorial
cd $SUBJECTS_DIR

setenv SUBJECTS_DIR $TUTORIAL_DATA/group_analysis_tutorial
cd $SUBJECTS_DIR

#qdec &

cd $SUBJECTS_DIR
mris_anatomical_stats -l lh.supramarg.label \
  -t lh.thickness -b -f 004/stats/lh.supramarg.stats 004 lh

#tkmedit 004 brainmask.mgz -label lh.supramarg.label

#------------------------------
# ---- AnatomicalROI ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA/group_analysis_tutorial
cd $SUBJECTS_DIR

#tkmedit 004 orig.mgz -aux aparc+aseg.mgz \
#  -seg aparc+aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt

#tksurfer 004 lh inflated -annot aparc.annot

cat $FREESURFER_HOME/FreeSurferColorLUT.txt

cat 004/label/lh.BA45.label

#tkmedit 004 orig.mgz

#tksurfer 004 lh inflated

cd $SUBJECTS_DIR/004/stats
cat aseg.stats

cd $SUBJECTS_DIR/004/stats
cat lh.aparc.stats

setenv SUBJECTS_DIR $TUTORIAL_DATA/group_analysis_tutorial
cd $SUBJECTS_DIR

asegstats2table --subjects 004 021 040 067 080 092 \
  --segno 11 17 18 \
  --tablefile aseg.vol.table

cat $FREESURFER_HOME/FreeSurferColorLUT.txt

cat aseg.vol.table

#oocalc aseg.vol.table

#/Applications/OpenOffice.org.app/Contents/MacOS/scalc aseg.vol.table

asegstats2table \
  --subjects 004 021 040 067 080 092 \
  --segno 11 17 18 \
  --meas mean \
  --tablefile aseg.mean-intensity.table

cat $FREESURFER_HOME/FreeSurferColorLUT.txt

asegstats2table \
  --subjects 004 021 040 067 080 092 \
  --segno 3007 3021 3022 4022 \
  --stats wmparc.stats \
  --tablefile wmparc.vol.table

cat $FREESURFER_HOME/FreeSurferColorLUT.txt

aparcstats2table --hemi lh \
  --subjects 004 021 040 067 080 092 \
  --tablefile lh.aparc.area.table

aparcstats2table --hemi lh \
  --subjects 004 021 040 067 080 092 \
  --meas thickness \
  --parc aparc.a2005s \
  --tablefile lh.aparc.a2005.thickness.table

#------------------------------
# ---- TroubleshootingData ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs

setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $SUBJECTS_DIR

#tkmedit topo_defect_before brainmask.mgz \
#  lh.white -aux wm.mgz -aux-surface rh.white
#
#
#tkmedit topo_defect_after brainmask.mgz \
#  lh.white -aux wm.mgz -aux-surface rh.white
#
#
#tkmedit wm1_edits_before brainmask.mgz \
#  lh.white -aux T1.mgz -aux-surface rh.white
#
#
#tksurfer wm1_edits_before lh inflated
#
#tkmedit pial_edits_before brainmask.mgz \
#  lh.white -aux T1.mgz -aux-surface rh.white
#
#
#tksurfer pial_edits_before lh inflated
#
#tkmedit skullstrip1_before brainmask.mgz \
#  lh.white -aux T1.mgz -aux-surface rh.white
#
#
#tkmedit cp_before brainmask.mgz \
#  lh.white -aux T1.mgz -aux-surface rh.white
#
#
#tksurfer cp_before lh inflated &
#tksurfer cp_before rh inflated &
#
#tkmedit tal_before brainmask.mgz \
#  lh.white -aux T1.mgz -aux-surface rh.white
#
#
#tksurfer tal_before lh inflated

#------------------------------
# ---- MultiModalRegistration ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $TUTORIAL_DATA/multimodal/fmri/fbirn-101

#tkregister2 --mov template.nii --s fbirn-anat-101.v4 \
# --regheader --reg myregister.dat --surf

#  bbregister --mov template.nii --bold \
#    --s fbirn-anat-101.v4 \
#    --init-fsl --reg register.dat

cat register.dat

#tkregister2 --mov template.nii --reg register.dat --surf

#------------------------------
# ---- MultiModalFmriIndividual ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $TUTORIAL_DATA/multimodal/fmri/fbirn-101

#tkregister2 --mov template.nii \
#  --reg bb.register.dat --surf

#tkmedit fbirn-anat-101.v4 orig.mgz -aux brain.mgz -seg aparc+aseg.mgz \
#  -overlay sig.nii -reg bb.register.dat -fthresh 2 -fmax 4

mri_vol2surf --mov sig.nii \
    --reg bb.register.dat \
    --projfrac 0.5 --interp nearest \
    --hemi lh --o lh.sig.mgh


mri_info lh.sig.mgh

#tksurfer fbirn-anat-101.v4 lh inflated -annot aparc \
#  -overlay lh.sig.mgh

mri_vol2vol --mov ces.nii \
    --reg bb.register.dat \
    --fstarg --interp nearest \
    --o ces.anat.bb.mgh

mri_info ces.anat.bb.mgh

mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.bb.stats

mri_vol2vol --mov sig.nii \
    --reg bb.register.dat \
    --fstarg --interp nearest \
    --o sig.anat.bb.mgh

mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.abs-masked.bb.stats \
   --mask sig.anat.bb.mgh --maskthresh 2 --masksign abs

mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.pos-masked.bb.stats \
   --mask sig.anat.bb.mgh --maskthresh 2 --masksign pos

#------------------------------
# ---- MultiModalFmriGroup ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $TUTORIAL_DATA/multimodal/fmri

#tkregister2 --mov fbirn-101/template.nii \
#  --reg fbirn-101/bb.register.dat --surf

mris_preproc --target fsaverage --hemi lh \
  --iv  fbirn-101/ces.nii fbirn-101/bb.register.dat \
  --iv  fbirn-103/ces.nii fbirn-103/bb.register.dat \
  --iv  fbirn-104/ces.nii fbirn-104/bb.register.dat \
  --iv  fbirn-105/ces.nii fbirn-105/bb.register.dat \
  --iv  fbirn-106/ces.nii fbirn-106/bb.register.dat \
  --projfrac 0.5 \
  --out lh.ces.mgh

mri_info lh.ces.mgh

mri_surf2surf --hemi lh --s fsaverage --fwhm 5 --cortex\
  --sval lh.ces.mgh --tval lh.ces.sm05.mgh

mri_glmfit --y lh.ces.sm05.mgh --surf fsaverage lh \
  --osgm --glmdir lh.ces.sm05.osgm --cortex

#tksurfer fsaverage lh inflated -annot aparc -ov lh.ces.sm05.osgm/osgm/sig.mgh

asegstats2table \
  --meas volume \
  --tablefile ces.pos-masked.vol.stats \
  --i fbirn-101/ces.pos-masked.bb.stats \
  fbirn-103/ces.pos-masked.bb.stats \
  fbirn-104/ces.pos-masked.bb.stats \
  fbirn-105/ces.pos-masked.bb.stats \
  fbirn-106/ces.pos-masked.bb.stats

asegstats2table \
  --meas mean \
  --tablefile ces.abs-masked.mean.stats \
  --i fbirn-101/ces.abs-masked.bb.stats \
  fbirn-103/ces.abs-masked.bb.stats \
  fbirn-104/ces.abs-masked.bb.stats \
  fbirn-105/ces.abs-masked.bb.stats \
  fbirn-106/ces.abs-masked.bb.stats

asegstats2table \
  --meas mean \
  --tablefile ces.pos-masked.mean.stats \
  --i fbirn-101/ces.pos-masked.bb.stats \
  fbirn-103/ces.pos-masked.bb.stats \
  fbirn-104/ces.pos-masked.bb.stats \
  fbirn-105/ces.pos-masked.bb.stats \
  fbirn-106/ces.pos-masked.bb.stats

#------------------------------
# ---- MultiModalDtiIndividual ----
#------------------------------

setenv TUTORIAL_DATA $FREESURFER_HOME/subjects/buckner_data/tutorial_subjs
setenv SUBJECTS_DIR $TUTORIAL_DATA
cd $TUTORIAL_DATA/multimodal/dti

#tkregister2 --mov lowb.nii --reg register.dat --surf

#tkmedit M87102113.v4 orig.mgz -aux brain.mgz \
#  -seg wmparc.mgz \
#  -reg register.dat -overlay fa.nii \
#  -fthresh 0.2 -fmax 1

mri_vol2vol --mov fa.nii \
  --reg register.dat \
  --fstarg --interp nearest \
  --o fa.anat.mgh

mri_segstats \
  --seg $SUBJECTS_DIR/M87102113.v4/mri/wmparc.mgz \
  --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
  --id 251 --id 3021 --id 3024 --id 3030 --id 12 --id 4 \
  --i fa.anat.mgh --sum fa.stats

#------------------------------
#Testing tksurfer, tkmedit and qdec
#------------------------------

cd $SUBJECTS_DIR

tkmedit good_output brainmask.mgz lh.white \
        -aux T1.mgz -aux-surface rh.white \
        -segmentation aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt

tksurfer good_output lh inflated

qdec 
