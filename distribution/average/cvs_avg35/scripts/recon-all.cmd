
 mri_convert /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun.mgz /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Tue Apr 20 09:52:46 EDT 2010

 cp /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/orig/001.mgz /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/rawavg.mgz 


 mri_convert /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/rawavg.mgz /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/transforms/talairach.xfm /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/orig.mgz /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/orig.mgz 

#--------------------------------------------
#@# Nu Intensity Correction Tue Apr 20 09:52:59 EDT 2010

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --n 2 

#--------------------------------------------
#@# Talairach Tue Apr 20 09:55:27 EDT 2010

 talairach_avi --i nu.mgz --xfm transforms/talairach.auto.xfm 


 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Tue Apr 20 09:56:19 EDT 2010

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /usr/local/freesurfer/dev/bin/extract_talairach_avi_QA.awk /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Intensity Normalization Tue Apr 20 09:56:20 EDT 2010

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Tue Apr 20 09:58:31 EDT 2010

 mri_em_register -skull nu.mgz /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 


 mri_watershed -T1 -brain_atlas /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 mri_gcut -110 -mult brainmask.auto.mgz T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Tue Apr 20 10:17:14 EDT 2010

 mri_em_register -mask brainmask.mgz nu.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Tue Apr 20 10:44:28 EDT 2010

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Tue Apr 20 10:46:20 EDT 2010

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Tue Apr 20 14:15:54 EDT 2010

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Tue Apr 20 14:16:55 EDT 2010

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Tue Apr 20 14:18:03 EDT 2010

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 

#--------------------------------------
#@# SubCort Seg Tue Apr 20 14:42:07 EDT 2010

 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/avg35rerun/mri/transforms/cc_up.lta avg35rerun 

#--------------------------------------
#@# Merge ASeg Tue Apr 20 15:02:48 EDT 2010

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Tue Apr 20 15:02:49 EDT 2010

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Tue Apr 20 15:06:20 EDT 2010

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Tue Apr 20 15:06:22 EDT 2010

 mri_segment brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Tue Apr 20 15:08:03 EDT 2010

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Tue Apr 20 15:08:59 EDT 2010

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Tue Apr 20 15:09:09 EDT 2010

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Tue Apr 20 15:09:14 EDT 2010

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Tue Apr 20 15:10:13 EDT 2010

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Tue Apr 20 15:17:06 EDT 2010

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 avg35rerun lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make Final Surf lh Tue Apr 20 15:28:40 EDT 2010

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs avg35rerun lh 

#--------------------------------------------
#@# Surf Volume lh Tue Apr 20 15:43:03 EDT 2010

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#--------------------------------------------
#@# Smooth2 lh Tue Apr 20 15:43:04 EDT 2010

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Tue Apr 20 15:43:09 EDT 2010

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Tue Apr 20 15:45:35 EDT 2010

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm avg35rerun lh curv sulc 

#--------------------------------------------
#@# Sphere lh Tue Apr 20 15:45:42 EDT 2010

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Tue Apr 20 17:12:59 EDT 2010

 mris_register -curv ../surf/lh.sphere /usr/local/freesurfer/dev/average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Tue Apr 20 18:04:26 EDT 2010

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Tue Apr 20 18:04:29 EDT 2010

 mrisp_paint -a 5 /usr/local/freesurfer/dev/average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Tue Apr 20 18:04:32 EDT 2010

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 avg35rerun lh ../surf/lh.sphere.reg /usr/local/freesurfer/dev/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats lh Tue Apr 20 18:05:26 EDT 2010

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab avg35rerun lh 

#-----------------------------------------
#@# Cortical Parc 2 lh Tue Apr 20 18:05:42 EDT 2010

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 avg35rerun lh ../surf/lh.sphere.reg /usr/local/freesurfer/dev/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Tue Apr 20 18:06:43 EDT 2010

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab avg35rerun lh 

#--------------------------------------------
#@# Tessellate rh Tue Apr 20 18:07:02 EDT 2010

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Tue Apr 20 18:07:11 EDT 2010

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Tue Apr 20 18:07:16 EDT 2010

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Tue Apr 20 18:08:13 EDT 2010

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Tue Apr 20 18:14:10 EDT 2010

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 avg35rerun rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make Final Surf rh Tue Apr 20 18:26:12 EDT 2010

 mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs avg35rerun rh 

#--------------------------------------------
#@# Surf Volume rh Tue Apr 20 18:40:24 EDT 2010

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#--------------------------------------------
#@# Smooth2 rh Tue Apr 20 18:40:25 EDT 2010

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Tue Apr 20 18:40:30 EDT 2010

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Tue Apr 20 18:42:51 EDT 2010

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm avg35rerun rh curv sulc 

#--------------------------------------------
#@# Sphere rh Tue Apr 20 18:42:57 EDT 2010

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Tue Apr 20 19:51:35 EDT 2010

 mris_register -curv ../surf/rh.sphere /usr/local/freesurfer/dev/average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Tue Apr 20 20:48:41 EDT 2010

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Tue Apr 20 20:48:43 EDT 2010

 mrisp_paint -a 5 /usr/local/freesurfer/dev/average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Tue Apr 20 20:48:46 EDT 2010

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 avg35rerun rh ../surf/rh.sphere.reg /usr/local/freesurfer/dev/average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/rh.aparc.annot 

#-----------------------------------------
#@# Parcellation Stats rh Tue Apr 20 20:49:39 EDT 2010

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab avg35rerun rh 

#-----------------------------------------
#@# Cortical Parc 2 rh Tue Apr 20 20:49:54 EDT 2010

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 avg35rerun rh ../surf/rh.sphere.reg /usr/local/freesurfer/dev/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Tue Apr 20 20:50:55 EDT 2010

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab avg35rerun rh 

#--------------------------------------------
#@# Cortical ribbon mask Tue Apr 20 20:51:12 EDT 2010

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon --save_distance avg35rerun 

#--------------------------------------------
#@# ASeg Stats Tue Apr 20 21:05:30 EDT 2010

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --ctab /usr/local/freesurfer/dev/ASegStatsLUT.txt --subject avg35rerun 

#-----------------------------------------
#@# AParc-to-ASeg Tue Apr 20 21:15:27 EDT 2010

 mri_aparc2aseg --s avg35rerun --volmask 


 mri_aparc2aseg --s avg35rerun --volmask --a2009s 

#-----------------------------------------
#@# WMParc Tue Apr 20 21:18:04 EDT 2010

 mri_aparc2aseg --s avg35rerun --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject avg35rerun --surf-wm-vol --ctab /usr/local/freesurfer/dev/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA Labels lh Tue Apr 20 21:36:58 EDT 2010
INFO: fsaverage subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to fsaverage subject...

 ln -s /usr/local/freesurfer/dev/subjects/fsaverage /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA1.label --trgsubject avg35rerun --trglabel ./lh.BA1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA2.label --trgsubject avg35rerun --trglabel ./lh.BA2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA3a.label --trgsubject avg35rerun --trglabel ./lh.BA3a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA3b.label --trgsubject avg35rerun --trglabel ./lh.BA3b.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA4a.label --trgsubject avg35rerun --trglabel ./lh.BA4a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA4p.label --trgsubject avg35rerun --trglabel ./lh.BA4p.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA6.label --trgsubject avg35rerun --trglabel ./lh.BA6.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA44.label --trgsubject avg35rerun --trglabel ./lh.BA44.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.BA45.label --trgsubject avg35rerun --trglabel ./lh.BA45.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.V1.label --trgsubject avg35rerun --trglabel ./lh.V1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.V2.label --trgsubject avg35rerun --trglabel ./lh.V2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/lh.MT.label --trgsubject avg35rerun --trglabel ./lh.MT.label --hemi lh --regmethod surface 


 mris_label2annot --s avg35rerun --hemi lh --ctab /usr/local/freesurfer/dev/average/colortable_BA.txt --l lh.BA1.label --l lh.BA2.label --l lh.BA3a.label --l lh.BA3b.label --l lh.BA4a.label --l lh.BA4p.label --l lh.BA6.label --l lh.BA44.label --l lh.BA45.label --l lh.V1.label --l lh.V2.label --l lh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/lh.BA.stats -b -a ./lh.BA.annot -c ./BA.ctab avg35rerun lh 

#--------------------------------------------
#@# BA Labels rh Tue Apr 20 21:39:09 EDT 2010

 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA1.label --trgsubject avg35rerun --trglabel ./rh.BA1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA2.label --trgsubject avg35rerun --trglabel ./rh.BA2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA3a.label --trgsubject avg35rerun --trglabel ./rh.BA3a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA3b.label --trgsubject avg35rerun --trglabel ./rh.BA3b.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA4a.label --trgsubject avg35rerun --trglabel ./rh.BA4a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA4p.label --trgsubject avg35rerun --trglabel ./rh.BA4p.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA6.label --trgsubject avg35rerun --trglabel ./rh.BA6.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA44.label --trgsubject avg35rerun --trglabel ./rh.BA44.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.BA45.label --trgsubject avg35rerun --trglabel ./rh.BA45.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.V1.label --trgsubject avg35rerun --trglabel ./rh.V1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.V2.label --trgsubject avg35rerun --trglabel ./rh.V2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/fsaverage/label/rh.MT.label --trgsubject avg35rerun --trglabel ./rh.MT.label --hemi rh --regmethod surface 


 mris_label2annot --s avg35rerun --hemi rh --ctab /usr/local/freesurfer/dev/average/colortable_BA.txt --l rh.BA1.label --l rh.BA2.label --l rh.BA3a.label --l rh.BA3b.label --l rh.BA4a.label --l rh.BA4p.label --l rh.BA6.label --l rh.BA44.label --l rh.BA45.label --l rh.V1.label --l rh.V2.label --l rh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/rh.BA.stats -b -a ./rh.BA.annot -c ./BA.ctab avg35rerun rh 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label lh Tue Apr 20 21:41:17 EDT 2010
INFO: lh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to lh.EC_average subject...

 ln -s /usr/local/freesurfer/dev/subjects/lh.EC_average /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/lh.EC_average 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o avg35rerun label lh.entorhinal lh sphere.reg lh.EC_average lh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/lh.entorhinal_exvivo.stats -b -l ./lh.entorhinal_exvivo.label avg35rerun lh 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label rh Tue Apr 20 21:41:35 EDT 2010
INFO: rh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to rh.EC_average subject...

 ln -s /usr/local/freesurfer/dev/subjects/rh.EC_average /autofs/cluster/con_001/users/lilla/CVS_atlas/buckner35rerun/rh.EC_average 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o avg35rerun label rh.entorhinal rh sphere.reg rh.EC_average rh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/rh.entorhinal_exvivo.stats -b -l ./rh.entorhinal_exvivo.label avg35rerun rh 

