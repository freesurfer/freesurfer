
 mri_convert /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152.mgz /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Fri Mar 30 18:25:37 EDT 2012

 cp /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/orig/001.mgz /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/rawavg.mgz 


 mri_convert /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/rawavg.mgz /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/orig.mgz -rt cubic --conform 


 mri_add_xform_to_header -c /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/transforms/talairach.xfm /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/orig.mgz /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Fri Mar 30 18:25:44 EDT 2012

 mri_nu_correct.mni --n 1 --proto-iters 1000 --distance 50 --no-rescale --i orig.mgz --o orig_nu.mgz 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 


 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Fri Mar 30 18:28:19 EDT 2012

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /usr/local/freesurfer/dev/bin/extract_talairach_avi_QA.awk /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Fri Mar 30 18:28:20 EDT 2012

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 

#--------------------------------------------
#@# Intensity Normalization Fri Mar 30 18:31:16 EDT 2012

 mri_normalize -g 1 nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Fri Mar 30 18:33:35 EDT 2012

 mri_em_register -skull nu.mgz /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta 


 mri_watershed -T1 -brain_atlas /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Fri Mar 30 19:03:03 EDT 2012

 mri_em_register -uns 3 -mask brainmask.mgz nu.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Fri Mar 30 19:19:27 EDT 2012

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Fri Mar 30 19:21:16 EDT 2012

 mri_ca_register -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca transforms/talairach.m3z 

#--------------------------------------
#@# CA Reg Inv Fri Mar 30 23:18:07 EDT 2012

 mri_ca_register -invert-and-save transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Fri Mar 30 23:19:30 EDT 2012

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca nu_noneck.mgz 

#--------------------------------------
#@# SkullLTA Fri Mar 30 23:20:42 EDT 2012

 mri_em_register -skull -t transforms/talairach.lta nu_noneck.mgz /usr/local/freesurfer/dev/average/RB_all_withskull_2008-03-26.gca transforms/talairach_with_skull_2.lta 

#--------------------------------------
#@# SubCort Seg Fri Mar 30 23:45:28 EDT 2012

 mri_ca_label -align -nobigventricles norm.mgz transforms/talairach.m3z /usr/local/freesurfer/dev/average/RB_all_2008-03-26.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/cvs_avg35_inMNI152/mri/transforms/cc_up.lta cvs_avg35_inMNI152 

#--------------------------------------
#@# Merge ASeg Sat Mar 31 00:06:43 EDT 2012

 cp aseg.auto.mgz aseg.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Sat Mar 31 00:06:44 EDT 2012

 mri_normalize -aseg aseg.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Sat Mar 31 00:10:29 EDT 2012

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Sat Mar 31 00:10:31 EDT 2012

 mri_segment brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Sat Mar 31 00:13:10 EDT 2012

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Sat Mar 31 00:14:05 EDT 2012

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Sat Mar 31 00:14:14 EDT 2012

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Sat Mar 31 00:14:20 EDT 2012

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Sat Mar 31 00:15:09 EDT 2012

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology lh Sat Mar 31 00:21:13 EDT 2012

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 cvs_avg35_inMNI152 lh 


 mris_euler_number ../surf/lh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 

#--------------------------------------------
#@# Make White Surf lh Sat Mar 31 00:36:57 EDT 2012

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs cvs_avg35_inMNI152 lh 

#--------------------------------------------
#@# Smooth2 lh Sat Mar 31 00:41:55 EDT 2012

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Sat Mar 31 00:42:00 EDT 2012

 mris_inflate ../surf/lh.smoothwm ../surf/lh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/lh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Sat Mar 31 00:44:23 EDT 2012

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm cvs_avg35_inMNI152 lh curv sulc 

#--------------------------------------------
#@# Sphere lh Sat Mar 31 00:44:31 EDT 2012

 mris_sphere -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Surf Reg lh Sat Mar 31 01:34:08 EDT 2012

 mris_register -curv ../surf/lh.sphere /usr/local/freesurfer/dev/average/lh.average.curvature.filled.buckner40.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Sat Mar 31 02:41:52 EDT 2012

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Sat Mar 31 02:41:55 EDT 2012

 mrisp_paint -a 5 /usr/local/freesurfer/dev/average/lh.average.curvature.filled.buckner40.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Sat Mar 31 02:41:57 EDT 2012

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 cvs_avg35_inMNI152 lh ../surf/lh.sphere.reg /usr/local/freesurfer/dev/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/lh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Sat Mar 31 02:42:50 EDT 2012

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs cvs_avg35_inMNI152 lh 

#--------------------------------------------
#@# Surf Volume lh Sat Mar 31 02:54:35 EDT 2012

 mris_calc -o lh.area.mid lh.area add lh.area.pial 


 mris_calc -o lh.area.mid lh.area.mid div 2 


 mris_calc -o lh.volume lh.area.mid mul lh.thickness 

#-----------------------------------------
#@# WM/GM Contrast lh Sat Mar 31 02:54:36 EDT 2012

 pctsurfcon --s cvs_avg35_inMNI152 --lh-only 

#-----------------------------------------
#@# Parcellation Stats lh Sat Mar 31 02:54:45 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab cvs_avg35_inMNI152 lh white 

#-----------------------------------------
#@# Cortical Parc 2 lh Sat Mar 31 02:55:01 EDT 2012

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 cvs_avg35_inMNI152 lh ../surf/lh.sphere.reg /usr/local/freesurfer/dev/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Sat Mar 31 02:55:59 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab cvs_avg35_inMNI152 lh white 

#--------------------------------------------
#@# Tessellate rh Sat Mar 31 02:56:17 EDT 2012

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 rh Sat Mar 31 02:56:27 EDT 2012

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 rh Sat Mar 31 02:56:32 EDT 2012

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere rh Sat Mar 31 02:57:18 EDT 2012

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology rh Sat Mar 31 03:03:42 EDT 2012

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 


 mris_fix_topology -mgz -sphere qsphere.nofix -ga -seed 1234 cvs_avg35_inMNI152 rh 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf rh Sat Mar 31 03:09:07 EDT 2012

 mris_make_surfaces -noaparc -whiteonly -mgz -T1 brain.finalsurfs cvs_avg35_inMNI152 rh 

#--------------------------------------------
#@# Smooth2 rh Sat Mar 31 03:14:05 EDT 2012

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 rh Sat Mar 31 03:14:10 EDT 2012

 mris_inflate ../surf/rh.smoothwm ../surf/rh.inflated 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 ../surf/rh.inflated 


#-----------------------------------------
#@# Curvature Stats rh Sat Mar 31 03:16:20 EDT 2012

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm cvs_avg35_inMNI152 rh curv sulc 

#--------------------------------------------
#@# Sphere rh Sat Mar 31 03:16:27 EDT 2012

 mris_sphere -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg rh Sat Mar 31 04:08:25 EDT 2012

 mris_register -curv ../surf/rh.sphere /usr/local/freesurfer/dev/average/rh.average.curvature.filled.buckner40.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white rh Sat Mar 31 05:16:58 EDT 2012

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv rh Sat Mar 31 05:17:00 EDT 2012

 mrisp_paint -a 5 /usr/local/freesurfer/dev/average/rh.average.curvature.filled.buckner40.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc rh Sat Mar 31 05:17:03 EDT 2012

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 cvs_avg35_inMNI152 rh ../surf/rh.sphere.reg /usr/local/freesurfer/dev/average/rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf rh Sat Mar 31 05:17:53 EDT 2012

 mris_make_surfaces -white NOWRITE -mgz -T1 brain.finalsurfs cvs_avg35_inMNI152 rh 

#--------------------------------------------
#@# Surf Volume rh Sat Mar 31 05:28:56 EDT 2012

 mris_calc -o rh.area.mid rh.area add rh.area.pial 


 mris_calc -o rh.area.mid rh.area.mid div 2 


 mris_calc -o rh.volume rh.area.mid mul rh.thickness 

#-----------------------------------------
#@# WM/GM Contrast rh Sat Mar 31 05:28:57 EDT 2012

 pctsurfcon --s cvs_avg35_inMNI152 --rh-only 

#-----------------------------------------
#@# Parcellation Stats rh Sat Mar 31 05:29:05 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab cvs_avg35_inMNI152 rh white 

#-----------------------------------------
#@# Cortical Parc 2 rh Sat Mar 31 05:29:21 EDT 2012

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 cvs_avg35_inMNI152 rh ../surf/rh.sphere.reg /usr/local/freesurfer/dev/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 rh Sat Mar 31 05:30:16 EDT 2012

 mris_anatomical_stats -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab cvs_avg35_inMNI152 rh white 

#--------------------------------------------
#@# Cortical ribbon mask Sat Mar 31 05:30:33 EDT 2012

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon cvs_avg35_inMNI152 

#--------------------------------------------
#@# ASeg Stats Sat Mar 31 05:41:25 EDT 2012

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --ctab /usr/local/freesurfer/dev/ASegStatsLUT.txt --subject cvs_avg35_inMNI152 

#-----------------------------------------
#@# AParc-to-ASeg Sat Mar 31 05:52:35 EDT 2012

 mri_aparc2aseg --s cvs_avg35_inMNI152 --volmask 


 mri_aparc2aseg --s cvs_avg35_inMNI152 --volmask --a2009s 

#-----------------------------------------
#@# WMParc Sat Mar 31 05:55:13 EDT 2012

 mri_aparc2aseg --s cvs_avg35_inMNI152 --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brain-vol-from-seg --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject cvs_avg35_inMNI152 --surf-wm-vol --ctab /usr/local/freesurfer/dev/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA Labels lh Sat Mar 31 06:16:30 EDT 2012
INFO: fsaverage subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to fsaverage subject...

 cd /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI; ln -s /usr/local/freesurfer/dev/subjects/fsaverage; cd - 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA1.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA2.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA3a.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA3a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA3b.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA3b.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA4a.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA4a.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA4p.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA4p.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA6.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA6.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA44.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA44.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.BA45.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.BA45.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.V1.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.V1.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.V2.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.V2.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/lh.MT.label --trgsubject cvs_avg35_inMNI152 --trglabel ./lh.MT.label --hemi lh --regmethod surface 


 mris_label2annot --s cvs_avg35_inMNI152 --hemi lh --ctab /usr/local/freesurfer/dev/average/colortable_BA.txt --l lh.BA1.label --l lh.BA2.label --l lh.BA3a.label --l lh.BA3b.label --l lh.BA4a.label --l lh.BA4p.label --l lh.BA6.label --l lh.BA44.label --l lh.BA45.label --l lh.V1.label --l lh.V2.label --l lh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/lh.BA.stats -b -a ./lh.BA.annot -c ./BA.ctab cvs_avg35_inMNI152 lh white 

#--------------------------------------------
#@# BA Labels rh Sat Mar 31 06:18:45 EDT 2012

 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA1.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA2.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA3a.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA3a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA3b.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA3b.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA4a.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA4a.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA4p.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA4p.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA6.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA6.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA44.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA44.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.BA45.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.BA45.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.V1.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.V1.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.V2.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.V2.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI/fsaverage/label/rh.MT.label --trgsubject cvs_avg35_inMNI152 --trglabel ./rh.MT.label --hemi rh --regmethod surface 


 mris_label2annot --s cvs_avg35_inMNI152 --hemi rh --ctab /usr/local/freesurfer/dev/average/colortable_BA.txt --l rh.BA1.label --l rh.BA2.label --l rh.BA3a.label --l rh.BA3b.label --l rh.BA4a.label --l rh.BA4p.label --l rh.BA6.label --l rh.BA44.label --l rh.BA45.label --l rh.V1.label --l rh.V2.label --l rh.MT.label --a BA --maxstatwinner --noverbose 


 mris_anatomical_stats -mgz -f ../stats/rh.BA.stats -b -a ./rh.BA.annot -c ./BA.ctab cvs_avg35_inMNI152 rh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label lh Sat Mar 31 06:20:56 EDT 2012
INFO: lh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to lh.EC_average subject...

 cd /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI; ln -s /usr/local/freesurfer/dev/subjects/lh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o cvs_avg35_inMNI152 label lh.entorhinal lh sphere.reg lh.EC_average lh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/lh.entorhinal_exvivo.stats -b -l ./lh.entorhinal_exvivo.label cvs_avg35_inMNI152 lh white 

#--------------------------------------------
#@# Ex-vivo Entorhinal Cortex Label rh Sat Mar 31 06:21:15 EDT 2012
INFO: rh.EC_average subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to rh.EC_average subject...

 cd /autofs/space/erdos_003/users/lzollei/my_stuff/research/tests/CVSvsMNI; ln -s /usr/local/freesurfer/dev/subjects/rh.EC_average; cd - 


 mris_spherical_average -erode 1 -orig white -t 0.4 -o cvs_avg35_inMNI152 label rh.entorhinal rh sphere.reg rh.EC_average rh.entorhinal_exvivo.label 


 mris_anatomical_stats -mgz -f ../stats/rh.entorhinal_exvivo.stats -b -l ./rh.entorhinal_exvivo.label cvs_avg35_inMNI152 rh white 

