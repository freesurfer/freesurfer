#!/bin/bash

HELP="

  Usage: test_tutorials.sh <options>

  Optional Arguments:

	-all			: Do all tutorials
	-quick			: Performs a quick subset of commands

	-auto_quit_freeview	: Automatically closes freeview after opening 
 	-skip_all_guis		: Skips all commands that open a GUI
	-skip_tk_guis		: Skips commands that open a tk GUI (tkmedit, tksurfer, etc)
	-skip_qdec_guis		: Skips commands that open qdec

	-individual_subject	: Do Interaction with Individual Subject Data tutorial
	-troubleshooting	: Do Troubleshooting tutorial
	-group_analysis		: Do Group Analysis tutorial
	-qdec			: Do QDEC tutorial
	-longitudinal		: Do Longitudinal tutorial
	-roi_analysis		: Do ROI Analysis tutorial
	-diffusion		: Do Diffusion tutorial
	-tracula		: Do TRACULA tutorial
	-fsfast			: Do FSFASt tutorial
	-multimodal		: DO Mulimodal tutorial

  Examples:

	Run all tutorial commands:

	  $> test_tutorials.sh -all

	Run a quick test to ensure basic functionaity:
 
   	  $> test_tutorials.sh -quick
	
	Run diffusion tutorial and multimodal tutorial:

	  $> test_tutorials.sh -diffusion -multimodal
	
	Run all tutorial commands and automatically close freeview after opening:

	  $> test_tutorials.sh -all _freeview

	Run all tutorial commands but skip all GUI calls:

	  $> test_tutorials.sh -all -skip_all_guis
"

do_Interaction_with_Individual_Subject_Data_Tutorial=0
do_Troubleshooting_Tutorial=0
do_Group_Analysis_Tutorial=0
do_QDEC_Tutorial=0
do_Logitudinal_Tutorial=0
do_ROI_Analysis_Tutorial=0
do_Diffusion_Processing_Tutorial=0
do_TRACULA_Tutorial=0
do_FSFAST_Tutorial=0
do_MultiModal_Integration=0
do_Quick_Test=0
auto_quit_freeview=0
skip_all_guis=0
skip_tk_guis=0
skip_qdec_guis=0

if [[ $# -lt 1 ]]; then
  echo ""
  echo "  ERROR: Requires at least 1 argument."
  echo ""
  echo "${HELP}"
  exit 1
fi

while [[ $# > 0 ]]; do
  key="$1"
  case $key in
    -all)
	do_Interaction_with_Individual_Subject_Data_Tutorial=1
	do_Troubleshooting_Tutorial=1
	do_Group_Analysis_Tutorial=1
	do_QDEC_Tutorial=1
	do_Logitudinal_Tutorial=1
	do_ROI_Analysis_Tutorial=1
	do_Diffusion_Processing_Tutorial=1
	do_TRACULA_Tutorial=1
	do_FSFAST_Tutorial=1
	do_MultiModal_Integration=1
    	;;
    -quick)
	do_Quick_Test=1
    	;;
    -auto_quit_freeview)
	auto_quit_freeview=1
        ;;
    -skip_all_guis)
	skip_all_guis=1
	;;
    -skip_tk_guis)
	skip_tk_guis=1
	;;
    -skip_qdec_guis)
	skip_qdec_guis=1
	;;
    -individual_subject)
	do_Interaction_with_Individual_Subject_Data_Tutorial=1
	;;
    -troubleshooting)
	do_Troubleshooting_Tutorial=1
	;;
    -group_analysis)
	do_Group_Analysis_Tutorial=1 
	;;
    -qdec)
	do_QDEC_Tutorial=1
    	;;
    -longitudinal)
	do_Logitudinal_Tutorial=1
	;;	
    -roi_analysis)
	do_ROI_Analysis_Tutorial=1
	;;
    -diffusion)
	do_Diffusion_Processing_Tutorial=1
	;;
    -tracula)
	do_TRACULA_Tutorial=1
	;;
    -fsfast)
     	do_FSFAST_Tutorial=1
 	;;
    -multimodal)
	do_MultiModal_Integration=1
	;;
    *)
    	echo
    	echo "${HELP}"
    	exit 1
    	;;
  esac
  shift
done


run_cmd () {
  cmd="${1}"

  ## Auto quit freeview?
  if [[ $auto_quit_freeview -eq 1 ]] && [[ $cmd == "freeview "* ]]; then
    cmd="${cmd} -quit"
  fi

  echo ""
  if [[ $skip_all_guis -eq 1 ]] && [[ $cmd == freeview* || $cmd == qdec* || $cmd == tkmedit* || $cmd == tkregis* || $cmd == tksurf* ]]; then
    echo "##################################################"
    echo "SKIPPING: ${cmd}"
    echo "##################################################"
    status=0
  elif [[ $skip_tk_guis -eq 1 ]] && [[ $cmd == tkmedit* || $cmd == tkregis* || $cmd == tksurf* ]]; then
    echo "##################################################"
    echo "SKIPPING: ${cmd}"
    echo "##################################################"
    status=0
  elif [[ $skip_qdec_guis -eq 1 ]] && [[ $cmd == "qdec"* ]]; then
    echo "##################################################"
    echo "SKIPPING: ${cmd}"
    echo "##################################################"
    status=0
  else
    echo "##################################################"
    echo "pwd=`pwd`"
    echo $cmd
    echo "##################################################"
    source $FREESURFER_HOME/sources.sh
    eval $cmd
    status=$?
  fi
  echo ""

  
  if [ $status -ne 0 ]; then
    echo "ERROR: CMD = $cmd"
    exit 1
  fi
}

if [ $do_Quick_Test -eq 1 ]; then

  # Interaction with Individual Subject
  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs
  run_cmd "cd $SUBJECTS_DIR"
  
  run_cmd "freeview -v \
  good_output/mri/T1.mgz \
  good_output/mri/wm.mgz \
  good_output/mri/brainmask.mgz \
  good_output/mri/aseg.mgz:colormap=lut:opacity=0.2 \
  -f good_output/surf/lh.white:edgecolor=blue \
  good_output/surf/lh.pial:edgecolor=red \
  good_output/surf/rh.white:edgecolor=blue \
  good_output/surf/rh.pial:edgecolor=red"


  # Group Analysis Tutorial
  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial
  
  run_cmd "cd $SUBJECTS_DIR/glm"

  run_cmd "mris_preproc --fsgd gender_age.fsgd \
  --cache-in thickness.fwhm10.fsaverage \
  --target fsaverage --hemi lh \
  --out lh.gender_age.thickness.10.mgh"

  # Test qdec 
  run_cmd "qdec &"

  # Diffusion 
  export SUBJECTS_DIR=$TUTORIAL_DATA/diffusion_recons

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/dmrirc.tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dwi_motion.txt"

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri"

  run_cmd "freeview $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dtifit_FA.nii.gz"

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/mni/"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.overall.txt"

  run_cmd "tractstats2table \
  --inputs $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.overall.txt \
  --overall --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.elmo.2012.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.elmo.2012.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list"

  run_cmd "tractstats2table --load-pathstats-from-file \
  $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list --overall \
  --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "tractstats2table --load-pathstats-from-file \
  $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list --overall \
  --only-measures FA_Avg --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.FA_Avg.table"
  
  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"

  run_cmd "trac-all -stat -c $TUTORIAL_DATA/diffusion_tutorial/dmrirc.tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"
 
  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"


  # FSFAST
  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/fsfast-tutorial.subjects"
  
  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"

  run_cmd "preproc-sess -s sess01 -fsd bold \
  -stc up -surface fsaverage lhrh -mni305 -fwhm 5 -per-run"

  # 
  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs"
  run_cmd "cd $SUBJECTS_DIR/multimodal/fmri"
  run_cmd "tkregister2 --mov fbirn-101/template.nii \
  --reg fbirn-101/bb.register.lta --surf &"

  run_cmd "mris_preproc --target fsaverage --hemi lh \
  --iv  fbirn-101/ces.nii fbirn-101/bb.register.lta \
  --iv  fbirn-103/ces.nii fbirn-103/bb.register.lta \
  --iv  fbirn-104/ces.nii fbirn-104/bb.register.lta \
  --iv  fbirn-105/ces.nii fbirn-105/bb.register.lta \
  --iv  fbirn-106/ces.nii fbirn-106/bb.register.lta \
  --projfrac 0.5 \
  --out lh.ces.mgh"

  run_cmd "mri_info lh.ces.mgh"

  run_cmd "mri_surf2surf --hemi lh --s fsaverage --fwhm 5 --cortex\
  --sval lh.ces.mgh --tval lh.ces.sm05.mgh"

  run_cmd "mri_glmfit --y lh.ces.sm05.mgh --surf fsaverage lh \
  --osgm --glmdir lh.ces.sm05.osgm --cortex"

  run_cmd "freeview -f \
  $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay=lh.ces.sm05.osgm/osgm/sig.mgh:overlay_threshold=2,5 \
  -viewport 3d"
 
  echo ""
  echo "  CONGRATULATIONS!! You passed the quick test."
  echo ""

fi

if [ $do_Interaction_with_Individual_Subject_Data_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs
  cd $SUBJECTS_DIR
  
  run_cmd "freeview -v \
  good_output/mri/T1.mgz \
  good_output/mri/wm.mgz \
  good_output/mri/brainmask.mgz \
  good_output/mri/aseg.mgz:colormap=lut:opacity=0.2 \
  -f good_output/surf/lh.white:edgecolor=blue \
  good_output/surf/lh.pial:edgecolor=red \
  good_output/surf/rh.white:edgecolor=blue \
  good_output/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -f  \
  good_output/surf/lh.pial:annot=aparc.annot:name=pial_aparc:visible=0 \
  good_output/surf/lh.inflated:overlay=lh.thickness:overlay_threshold=0.1,3::name=inflated_thickness:visible=0 \
  good_output/surf/lh.inflated:visible=0 \
  good_output/surf/lh.white:visible=0 \
  good_output/surf/lh.pial \
  --viewport 3d"

  run_cmd "freeview -v \
  good_output/mri/T1.mgz \
  good_output/mri/wm.mgz \
  good_output/mri/brainmask.mgz \
  good_output/mri/aseg.mgz:colormap=lut:opacity=0.2 \
  -f good_output/surf/lh.white:edgecolor=blue \
  good_output/surf/lh.pial:edgecolor=red \
  good_output/surf/rh.white:edgecolor=blue \
  good_output/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -f good_output/surf/lh.pial:annot=aparc.annot:name=pial_aparc:visible=0 \
  good_output/surf/lh.inflated:overlay=lh.thickness:overlay_threshold=0.1,3::name=inflated_thickness:visible=0 \
  good_output/surf/lh.inflated:visible=0 \
  good_output/surf/lh.white:visible=0 \
  good_output/surf/lh.pial \
  --viewport 3d"

fi

if [ $do_Troubleshooting_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs
  
  run_cmd "cd $SUBJECTS_DIR"
  
  run_cmd "freeview -v  pial_edits_before/mri/T1.mgz  \
  pial_edits_before/mri/brainmask.mgz  \
  -f pial_edits_before/surf/lh.white:edgecolor=yellow \
  pial_edits_before/surf/lh.pial:edgecolor=red \
  pial_edits_before/surf/rh.white:edgecolor=yellow \
  pial_edits_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v pial_edits_before/mri/T1.mgz \
  pial_edits_before/mri/brainmask.mgz \
  -f pial_edits_before/surf/lh.white:edgecolor=yellow \
  pial_edits_before/surf/lh.pial:edgecolor=red \
  pial_edits_before/surf/rh.white:edgecolor=yellow \
  pial_edits_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v wm1_edits_before/mri/brainmask.mgz \
  wm1_edits_before/mri/wm.mgz:colormap=heat:opacity=0.4 \
  -f wm1_edits_before/surf/lh.white:edgecolor=blue \
  wm1_edits_before/surf/lh.pial:edgecolor=red \
  wm1_edits_before/surf/rh.white:edgecolor=blue \
  wm1_edits_before/surf/rh.pial:edgecolor=red \
  wm1_edits_before/surf/rh.inflated:visible=0 \
  wm1_edits_before/surf/lh.inflated:visible=0"

  run_cmd "freeview -v wm1_edits_before/mri/brainmask.mgz \
  wm1_edits_before/mri/wm.mgz:colormap=heat:opacity=0.4 \
  -f wm1_edits_before/surf/lh.white:edgecolor=blue \
  wm1_edits_before/surf/lh.pial:edgecolor=red \
  wm1_edits_before/surf/rh.white:edgecolor=blue \
  wm1_edits_before/surf/rh.pial:edgecolor=red \
  wm1_edits_before/surf/rh.inflated:visible=0 \
  wm1_edits_before/surf/lh.inflated:visible=0"

  run_cmd "freeview -v wm1_edits_after/mri/brainmask.mgz \
  wm1_edits_after/mri/wm.mgz:colormap=heat:opacity=0.4 \
  -f wm1_edits_after/surf/lh.white:edgecolor=blue \
  wm1_edits_after/surf/lh.pial:edgecolor=red \
  wm1_edits_after/surf/rh.white:edgecolor=blue \
  wm1_edits_after/surf/rh.pial:edgecolor=red \
  wm1_edits_after/surf/rh.inflated:visible=0 \
  wm1_edits_after/surf/lh.inflated:visible=0"

  run_cmd "freeview -v wm1_edits_after/mri/T1.mgz  \
  wm1_edits_after/mri/brainmask.mgz  \
  -f wm1_edits_after/surf/lh.white:edgecolor=yellow \
  wm1_edits_after/surf/lh.pial:edgecolor=red \
  wm1_edits_after/surf/rh.white:edgecolor=yellow \
  wm1_edits_after/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v topo_defect_before/mri/brainmask.mgz \
  topo_defect_before/mri/wm.mgz:colormap=heat:opacity=0.4 \
  -f topo_defect_before/surf/lh.white:edgecolor=yellow \
  topo_defect_before/surf/lh.pial:edgecolor=red \
  topo_defect_before/surf/rh.white:edgecolor=yellow \
  topo_defect_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v topo_defect_before/mri/brainmask.mgz \
  topo_defect_before/mri/wm.mgz:colormap=heat:opacity=0.4:visible=0 \
  -f topo_defect_before/surf/lh.white:edgecolor=yellow \
  topo_defect_before/surf/lh.pial:edgecolor=red \
  topo_defect_before/surf/rh.white:edgecolor=yellow \
  topo_defect_before/surf/rh.pial:edgecolor=red \
  topo_defect_before/surf/lh.smoothwm.nofix:visible=0"

  run_cmd "freeview -v topo_defect_after/mri/brainmask.mgz \
  topo_defect_after/mri/wm.mgz:colormap=heat:opacity=0.4:visible=0 \
  -f topo_defect_after/surf/lh.white:edgecolor=yellow \
  topo_defect_after/surf/lh.pial:edgecolor=red \
  topo_defect_after/surf/rh.white:edgecolor=yellow \
  topo_defect_after/surf/rh.pial:edgecolor=red \
  topo_defect_after/surf/lh.smoothwm.nofix:visible=0"

  run_cmd "freeview -v topo_defect_after/mri/brainmask.mgz \
  topo_defect_after/mri/wm.mgz:colormap=heat:opacity=0.4 \
  -f topo_defect_after/surf/lh.white:edgecolor=yellow \
  topo_defect_after/surf/lh.pial:edgecolor=red \
  topo_defect_after/surf/rh.white:edgecolor=yellow \
  topo_defect_after/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v skullstrip1_before/mri/T1.mgz \
  skullstrip1_before/mri/brainmask.mgz \
  -f skullstrip1_before/surf/lh.white:edgecolor=yellow \
  skullstrip1_before/surf/lh.pial:edgecolor=red \
  skullstrip1_before/surf/rh.white:edgecolor=yellow \
  skullstrip1_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v skullstrip1_before/mri/T1.mgz \
  skullstrip1_before/mri/brainmask.mgz \
  -f skullstrip1_before/surf/lh.white:edgecolor=yellow \
  skullstrip1_before/surf/lh.pial:edgecolor=red \
  skullstrip1_before/surf/rh.white:edgecolor=yellow \
  skullstrip1_before/surf/rh.pial:edgecolor=red"

  run_cmd "recon-all -skullstrip -wsthresh 35 -clean-bm -no-wsgcaatlas -subjid skullstrip1_before"

  run_cmd "freeview -v skullstrip1_after/mri/T1.mgz \
  skullstrip1_after/mri/brainmask.mgz \
  -f skullstrip1_after/surf/lh.white:edgecolor=yellow \
  skullstrip1_after/surf/lh.pial:edgecolor=red \
  skullstrip1_after/surf/rh.white:edgecolor=yellow \
  skullstrip1_after/surf/rh.pial:edgecolor=red"

  #run_cmd "cd $SUBJECTS_DIR/yoursubj/mri/transforms"
  ##run_cmd "cp talairach_with_skull.lta bak"
  #run_cmd "cp talairach_with_skull_2.lta talairach_with_skull.lta"
  #run_cmd "recon-all -s yoursubj -skullstrip -clean-bm -clean-lta"

  run_cmd "freeview -v skullstrip1_after/mri/T1.mgz \
  skullstrip1_after/mri/brainmask.mgz \
  -f skullstrip1_after/surf/lh.white:edgecolor=yellow \
  skullstrip1_after/surf/lh.pial:edgecolor=red \
  skullstrip1_after/surf/rh.white:edgecolor=yellow \
  skullstrip1_after/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v cp_before/mri/brainmask.mgz \
  cp_before/mri/T1.mgz \
  -f cp_before/surf/lh.white:edgecolor=blue \
  cp_before/surf/lh.pial:edgecolor=red \
  cp_before/surf/rh.white:edgecolor=blue \
  cp_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -v cp_before/mri/brainmask.mgz \
  -f cp_before/surf/lh.white:edgecolor=blue \
  cp_before/surf/lh.pial:edgecolor=red \
  cp_before/surf/rh.white:edgecolor=blue \
  cp_before/surf/rh.pial:edgecolor=red"

  run_cmd "freeview -f \
  tal_before/surf/lh.inflated:visible=0 \
  tal_before/surf/rh.inflated \
  -viewport 3d"

  run_cmd "freeview -v tal_before/mri/T1.mgz \
  -v tal_before/mri/brainmask.mgz:reg=tal_before/mri/transforms/talairach.xfm"

  run_cmd "tkregister2 --mgz --s tal_before --fstal --surf orig"

  run_cmd "freeview -v tal_before/mri/T1.mgz \
  tal_before/mri/brainmask.mgz:reg=tal_before/mri/transforms/talairach.xfm \
  tal_after/mri/brainmask.mgz:reg=tal_after/mri/transforms/talairach.xfm"

 # run_cmd "recon-all -s tal_before -talairach -use-mritotal -tal-check -clean-tal"

fi

if [ $do_Group_Analysis_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial
  
  run_cmd "cd $SUBJECTS_DIR/glm"

  run_cmd "mris_preproc --fsgd gender_age.fsgd \
  --cache-in thickness.fwhm10.fsaverage \
  --target fsaverage --hemi lh \
  --out lh.gender_age.thickness.10.mgh"

  run_cmd "mri_info lh.gender_age.thickness.10.mgh"

  ## Do Not Run
  #run_cmd "mris_preproc --fsgd gender_age.fsgd \
  #--target fsaverage --hemi lh \
  #--meas thickness \
  #--out lh.gender_age.thickness.00.mgh"

  run_cmd "mri_surf2surf --hemi lh \
  --s fsaverage \
  --sval lh.gender_age.thickness.00.mgh \
  --fwhm 10 \
  --cortex \
  --tval lh.gender_age.thickness.10B.mgh"

  run_cmd "mri_glmfit \
  --y lh.gender_age.thickness.10.mgh \
  --fsgd gender_age.fsgd dods\
  --C lh-Avg-thickness-age-Cor.mtx \
  --surf fsaverage lh \
  --cortex \
  --glmdir lh.gender_age.glmdir"

  run_cmd "freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay=lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/sig.mgh:overlay_threshold=4,5 -viewport 3d"

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial
  
  run_cmd "cd $SUBJECTS_DIR/glm"

  run_cmd "mri_glmfit-sim \
  --glmdir lh.gender_age.glmdir \
  --cache 4 neg \
  --cwp  0.05\
  --2spaces"

  run_cmd "cat lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/cache.th40.neg.sig.cluster.summary"

  run_cmd "freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:overlay=lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/cache.th40.neg.sig.cluster.mgh:overlay_threshold=2,5:annot=lh.gender_age.glmdir/lh-Avg-thickness-age-Cor/cache.th40.neg.sig.ocn.annot -viewport 3d"

fi

if [ $do_QDEC_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs
 
  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "if [ ! -e fsaverage ];then ln -s $FREESURFER_HOME/subjects/fsaverage;fi"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial"

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "qdec"

  run_cmd "mri_label2label --srclabel lh.supramarg --srcsubject fsaverage \
  --trgsubject 004 --trglabel lh.supramarg --regmethod surface --hemi lh"

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "mris_anatomical_stats -l lh.supramarg.label \
  -t lh.thickness -b -f 004/stats/lh.supramarg.stats 004 lh"

fi

if [ $do_Logitudinal_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/long-tutorial

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "freeview -v OAS2_0001/mri/norm.mgz \
  -f OAS2_0001/surf/lh.pial:edgecolor=red \
  OAS2_0001/surf/rh.pial:edgecolor=red \
  OAS2_0001/surf/lh.white:edgecolor=blue \
  OAS2_0001/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0001_MR1.long.OAS2_0001/mri/norm.mgz \
  OAS2_0001_MR2.long.OAS2_0001/mri/norm.mgz \
  -f OAS2_0001_MR1.long.OAS2_0001/surf/lh.pial:edgecolor=red \
  OAS2_0001_MR1.long.OAS2_0001/surf/lh.white:edgecolor=blue \
  OAS2_0001_MR2.long.OAS2_0001/surf/lh.pial:edgecolor=255,128,128 \
  OAS2_0001_MR2.long.OAS2_0001/surf/lh.white:edgecolor=lightblue"

  run_cmd "cat qdec/long.qdec.table.dat"

  ## DO NOT RUN THIS COMMANS. It will take a while and has already been done for you. 
  #run_cmd "long_mris_slopes --qdec ./qdec/long.qdec.table.dat \
  #--meas thickness --hemi lh --do-avg --do-rate --do-pc1 --do-spc \
  #--do-stack --do-label --time years --qcache fsaverage --sd $SUBJECTS_DIR"

  run_cmd "freeview -f \
  OAS2_0001/surf/lh.pial:overlay=OAS2_0001/surf/lh.long.thickness-avg.fwhm15.mgh:overlay_threshold=0,3.5:overlay=OAS2_0001/surf/lh.long.thickness-stack.mgh:annot=OAS2_0001/label/lh.aparc.annot:annot_outline=1 --timecourse --colorscale"

  run_cmd "freeview -f \
  fsaverage/surf/lh.pial:overlay=$SUBJECTS_DIR/OAS2_0001/surf/lh.long.thickness-spc.fwhm15.fsaverage.mgh:overlay_threshold=2,5"

  run_cmd "long_qdec_table --qdec ./qdec/long.qdec.table.dat \
  --cross --out ./qdec/cross.qdec.table.dat"

  run_cmd "qdec --table ./qdec/cross.qdec.table.dat"

  export SUBJECTS_DIR=$TUTORIAL_DATA/long-tutorial

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "freeview -v OAS2_0004/mri/T1.mgz \
  OAS2_0004/mri/brainmask.mgz \
  -f OAS2_0004/surf/lh.pial:edgecolor=red \
  OAS2_0004/surf/rh.pial:edgecolor=red \
  OAS2_0004/surf/lh.white:edgecolor=blue \
  OAS2_0004/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR1/mri/T1.mgz \
  OAS2_0004_MR1/mri/brainmask.mgz \
  -f OAS2_0004_MR1/surf/lh.pial:edgecolor=red \
  OAS2_0004_MR1/surf/rh.pial:edgecolor=red \
  OAS2_0004_MR1/surf/lh.white:edgecolor=blue \
  OAS2_0004_MR1/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR2/mri/T1.mgz \
  OAS2_0004_MR2/mri/brainmask.mgz \
  -f OAS2_0004_MR2/surf/lh.pial:edgecolor=red \
  OAS2_0004_MR2/surf/rh.pial:edgecolor=red \
  OAS2_0004_MR2/surf/lh.white:edgecolor=blue \
  OAS2_0004_MR2/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR1.long.OAS2_0004/mri/T1.mgz \
  OAS2_0004_MR1.long.OAS2_0004/mri/brainmask.mgz \
  -f OAS2_0004_MR1.long.OAS2_0004/surf/lh.pial:edgecolor=red \
  OAS2_0004_MR1.long.OAS2_0004/surf/rh.pial:edgecolor=red \
  OAS2_0004_MR1.long.OAS2_0004/surf/lh.white:edgecolor=blue \
  OAS2_0004_MR1.long.OAS2_0004/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR2.long.OAS2_0004/mri/T1.mgz \
  OAS2_0004_MR2.long.OAS2_0004/mri/brainmask.mgz \
  -f OAS2_0004_MR2.long.OAS2_0004/surf/lh.pial:edgecolor=red \
  OAS2_0004_MR2.long.OAS2_0004/surf/rh.pial:edgecolor=red \
  OAS2_0004_MR2.long.OAS2_0004/surf/lh.white:edgecolor=blue \
  OAS2_0004_MR2.long.OAS2_0004/surf/rh.white:edgecolor=blue"

  run_cmd "recon-all -subjid OAS2_0004_MR1 -skullstrip -wsthresh 20 -clean-bm -no-wsgcaatlas"

  run_cmd "recon-all -subjid OAS2_0004_MR2 -skullstrip -wsthresh 55 -clean-bm -no-wsgcaatlas"

 # run_cmd "recon-all -subjid OAS2_0004_MR1 -autorecon2 -autorecon3"

  #run_cmd "recon-all -subjid OAS2_0004_MR2 -autorecon2 -autorecon3"

  run_cmd "freeview -v OAS2_0004_MR1_fixed/mri/T1.mgz \
            OAS2_0004_MR1_fixed/mri/brainmask.mgz \
         -f OAS2_0004_MR1_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0004_MR1_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0004_MR1_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0004_MR1_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR2_fixed/mri/T1.mgz \
            OAS2_0004_MR2_fixed/mri/brainmask.mgz \
         -f OAS2_0004_MR2_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0004_MR2_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0004_MR2_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0004_MR2_fixed/surf/rh.white:edgecolor=blue"


  run_cmd "freeview -v OAS2_0057/mri/T1.mgz \
            OAS2_0057/mri/brainmask.mgz \
         -f OAS2_0057/surf/lh.pial:edgecolor=red \
            OAS2_0057/surf/rh.pial:edgecolor=red \
            OAS2_0057/surf/lh.white:edgecolor=blue \
            OAS2_0057/surf/rh.white:edgecolor=blue"

  #run_cmd "recon-all -skullstrip -subjid OAS2_0004_MR2"

  #run_cmd "recon-all -autorecon2 -autorecon3 -subjid OAS2_0004_MR2"

  #run_cmd "recon-all -long OAS2_0004_MR1 OAS2_0004 -all"

  #run_cmd "recon-all -long OAS2_0004_MR2 OAS2_0004 -all"

  run_cmd "freeview -v OAS2_0004_fixed/mri/T1.mgz \
            OAS2_0004_fixed/mri/brainmask.mgz \
         -f OAS2_0004_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0004_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0004_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0004_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR1.long.OAS2_0004_fixed/mri/T1.mgz \
            OAS2_0004_MR1.long.OAS2_0004_fixed/mri/brainmask.mgz \
         -f OAS2_0004_MR1.long.OAS2_0004_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0004_MR1.long.OAS2_0004_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0004_MR1.long.OAS2_0004_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0004_MR1.long.OAS2_0004_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0004_MR2.long.OAS2_0004_fixed/mri/T1.mgz \
            OAS2_0004_MR2.long.OAS2_0004_fixed/mri/brainmask.mgz \
         -f OAS2_0004_MR2.long.OAS2_0004_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0004_MR2.long.OAS2_0004_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0004_MR2.long.OAS2_0004_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0004_MR2.long.OAS2_0004_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0057_MR1/mri/T1.mgz \
            OAS2_0057_MR1/mri/brainmask.mgz \
         -f OAS2_0057_MR1/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR1/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR1/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR1/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0057_MR2/mri/T1.mgz \
            OAS2_0057_MR2/mri/brainmask.mgz \
         -f OAS2_0057_MR2/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR2/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR2/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR2/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0057_MR1.long.OAS2_0057/mri/T1.mgz \
            OAS2_0057_MR1.long.OAS2_0057/mri/brainmask.mgz \
         -f OAS2_0057_MR1.long.OAS2_0057/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR1.long.OAS2_0057/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0057_MR2.long.OAS2_0057/mri/T1.mgz \
            OAS2_0057_MR2.long.OAS2_0057/mri/brainmask.mgz \
         -f OAS2_0057_MR2.long.OAS2_0057/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR2.long.OAS2_0057/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR2.long.OAS2_0057/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR2.long.OAS2_0057/surf/rh.white:edgecolor=blue"

  #run_cmd "recon-all -base OAS2_0057 -autorecon-pial"

  #run_cmd "recon-all -long OAS2_0057_MR1 OAS2_0057 -all"

  #run_cmd "recon-all -long OAS2_0057_MR2 OAS2_0057 -all"

  run_cmd "freeview -v OAS2_0057_MR1.long.OAS2_0057/mri/T1.mgz \
            OAS2_0057_MR1.long.OAS2_0057/mri/brainmask.mgz \
         -f OAS2_0057_MR1.long.OAS2_0057/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR1.long.OAS2_0057/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0057_MR1.long.OAS2_0057_fixed/mri/T1.mgz \
            OAS2_0057_MR1.long.OAS2_0057_fixed/mri/brainmask.mgz \
         -f OAS2_0057_MR1.long.OAS2_0057_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0057_MR1.long.OAS2_0057_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0057_MR1.long.OAS2_0057_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR1/mri/T1.mgz \
            OAS2_0121_MR1/mri/brainmask.mgz \
         -f OAS2_0121_MR1/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR2/mri/T1.mgz \
            OAS2_0121_MR2/mri/brainmask.mgz \
         -f OAS2_0121_MR2/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR2/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR2/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR2/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121/mri/T1.mgz \
            OAS2_0121/mri/brainmask.mgz \
         -f OAS2_0121/surf/lh.pial:edgecolor=red \
            OAS2_0121/surf/rh.pial:edgecolor=red \
            OAS2_0121/surf/lh.white:edgecolor=blue \
            OAS2_0121/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR1.long.OAS2_0121/mri/T1.mgz \
            OAS2_0121_MR1.long.OAS2_0121/mri/brainmask.mgz \
         -f OAS2_0121_MR1.long.OAS2_0121/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1.long.OAS2_0121/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR2.long.OAS2_0121/mri/T1.mgz \
            OAS2_0121_MR2.long.OAS2_0121/mri/brainmask.mgz \
         -f OAS2_0121_MR2.long.OAS2_0121/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR2.long.OAS2_0121/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR2.long.OAS2_0121/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR2.long.OAS2_0121/surf/rh.white:edgecolor=blue"

 # run_cmd "recon-all -subjid OAS2_0121_MR1 -autorecon2 -autorecon3"

 # run_cmd "recon-all -subjid OAS2_0121_MR2 -autorecon2 -autorecon3"

  run_cmd "freeview -v OAS2_0121_MR1_fixed/mri/T1.mgz \
            OAS2_0121_MR1_fixed/mri/brainmask.mgz \
         -f OAS2_0121_MR1_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR2_fixed/mri/T1.mgz \
            OAS2_0121_MR2_fixed/mri/brainmask.mgz \
         -f OAS2_0121_MR2_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR2_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR2_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR2_fixed/surf/rh.white:edgecolor=blue"

 # run_cmd "recon-all -base OAS2_0121 -tp OAS2_0121_MR1 -tp OAS2_0121_MR2 -all"

  run_cmd "freeview -v OAS2_0121_intermediate/mri/T1.mgz \
            OAS2_0121_intermediate/mri/brainmask.mgz \
         -f OAS2_0121_intermediate/surf/lh.pial:edgecolor=red \
            OAS2_0121_intermediate/surf/rh.pial:edgecolor=red \
            OAS2_0121_intermediate/surf/lh.white:edgecolor=blue \
            OAS2_0121_intermediate/surf/rh.white:edgecolor=blue"

 # run_cmd "recon-all -base OAS2_0121 -autorecon2 -autorecon3"

  #run_cmd "recon-all -long OAS2_0121_MR1 OAS2_0121 -all"

 # run_cmd "recon-all -long OAS2_0121_MR2 OAS2_0121 -all"

  run_cmd "freeview -v OAS2_0121_MR1.long.OAS2_0121/mri/T1.mgz \
            OAS2_0121_MR1.long.OAS2_0121/mri/brainmask.mgz \
         -f OAS2_0121_MR1.long.OAS2_0121/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1.long.OAS2_0121/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR1.long.OAS2_0121_intermediate/mri/T1.mgz \
            OAS2_0121_MR1.long.OAS2_0121_intermediate/mri/brainmask.mgz \
         -f OAS2_0121_MR1.long.OAS2_0121_intermediate/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121_intermediate/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121_intermediate/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1.long.OAS2_0121_intermediate/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0121_MR1.long.OAS2_0121_fixed/mri/T1.mgz \
            OAS2_0121_MR1.long.OAS2_0121_fixed/mri/brainmask.mgz \
         -f OAS2_0121_MR1.long.OAS2_0121_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0121_MR1.long.OAS2_0121_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0121_MR1.long.OAS2_0121_fixed/surf/rh.white:edgecolor=blue"



  run_cmd "freeview -v OAS2_0185/mri/wm.mgz \
            OAS2_0185/mri/brainmask.mgz \
         -f OAS2_0185/surf/lh.pial:edgecolor=red \
            OAS2_0185/surf/rh.pial:edgecolor=red \
            OAS2_0185/surf/lh.white:edgecolor=blue \
            OAS2_0185/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0185_MR1/mri/wm.mgz \
            OAS2_0185_MR1/mri/brainmask.mgz \
         -f OAS2_0185_MR1/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR1/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR1/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR1/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0185_MR2/mri/wm.mgz \
            OAS2_0185_MR2/mri/brainmask.mgz \
         -f OAS2_0185_MR2/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR2/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR2/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR2/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0185_MR1.long.OAS2_0185/mri/wm.mgz \
            OAS2_0185_MR1.long.OAS2_0185/mri/brainmask.mgz \
         -f OAS2_0185_MR1.long.OAS2_0185/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR1.long.OAS2_0185/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR1.long.OAS2_0185/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR1.long.OAS2_0185/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0185_MR2.long.OAS2_0185/mri/wm.mgz \
            OAS2_0185_MR2.long.OAS2_0185/mri/brainmask.mgz \
         -f OAS2_0185_MR2.long.OAS2_0185/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR2.long.OAS2_0185/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR2.long.OAS2_0185/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR2.long.OAS2_0185/surf/rh.white:edgecolor=blue"

  #run_cmd "recon-all -base OAS2_0185 -autorecon2-wm -autorecon3"


  run_cmd "freeview -v OAS2_0185_fixed/mri/wm.mgz \
            OAS2_0185_fixed/mri/brainmask.mgz \
         -f OAS2_0185_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0185_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0185_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0185_fixed/surf/rh.white:edgecolor=blue"

  #run_cmd "recon-all -long OAS2_0185_MR1 OAS2_0185 -all"

#  run_cmd "recon-all -long OAS2_0185_MR2 OAS2_0185 -all"

  run_cmd "freeview -v OAS2_0185_MR1.long.OAS2_0185_fixed/mri/wm.mgz \
            OAS2_0185_MR1.long.OAS2_0185_fixed/mri/brainmask.mgz \
         -f OAS2_0185_MR1.long.OAS2_0185_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR1.long.OAS2_0185_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR1.long.OAS2_0185_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR1.long.OAS2_0185_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0185_MR2.long.OAS2_0185_fixed/mri/wm.mgz \
            OAS2_0185_MR2.long.OAS2_0185_fixed/mri/brainmask.mgz \
         -f OAS2_0185_MR2.long.OAS2_0185_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0185_MR2.long.OAS2_0185_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0185_MR2.long.OAS2_0185_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0185_MR2.long.OAS2_0185_fixed/surf/rh.white:edgecolor=blue"


  run_cmd "freeview -v OAS2_0002_MR2/mri/brain.finalsurfs.mgz \
            OAS2_0002_MR2/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR2/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR2/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR2/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR2/surf/rh.white:edgecolor=blue"

  run_cmd "cp OAS2_0002_MR2/mri/brain.finalsurfs.mgz OAS2_0002_MR2/mri/brain.finalsurfs.manedit.mgz"

#  run_cmd "recon-all -subjid OAS2_0002_MR2 -autorecon-pial"
  
  run_cmd "freeview -v OAS2_0002_MR2_fixed/mri/brainmask.mgz \
            OAS2_0002_MR2_fixed/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_MR2_fixed/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR2_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR2_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR2_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR2_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "cp OAS2_0002_MR1/mri/brain.finalsurfs.mgz OAS2_0002_MR1/mri/brain.finalsurfs.manedit.mgz"


  run_cmd "freeview -v OAS2_0002_MR1/mri/brainmask.mgz \
            OAS2_0002_MR1/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_MR1/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR1/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR1/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR1/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR1/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0002_MR1_fixed/mri/brainmask.mgz \
            OAS2_0002_MR1_fixed/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_MR1_fixed/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR1_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR1_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR1_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR1_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0002/mri/brain.finalsurfs.mgz \
            OAS2_0002/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002/surf/lh.pial:edgecolor=red \
            OAS2_0002/surf/rh.pial:edgecolor=red \
            OAS2_0002/surf/lh.white:edgecolor=blue \
            OAS2_0002/surf/rh.white:edgecolor=blue"

  run_cmd "cp OAS2_0002/mri/brain.finalsurfs.mgz OAS2_0002/mri/brain.finalsurfs.manedit.mgz"

 # run_cmd "recon-all -base OAS2_0002 -autorecon-pial"

  run_cmd "freeview -v OAS2_0002_fixed/mri/brainmask.mgz \
            OAS2_0002_fixed/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_fixed/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0002_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0002_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0002_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0002_MR2.long.OAS2_0002/mri/brain.finalsurfs.mgz \
            OAS2_0002_MR2.long.OAS2_0002/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR2.long.OAS2_0002/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR2.long.OAS2_0002/surf/rh.white:edgecolor=blue"

  run_cmd "cp OAS2_0002_MR2.long.OAS2_0002/mri/brain.finalsurfs.mgz OAS2_0002_MR2.long.OAS2_0002/mri/brain.finalsurfs.manedit.mgz"

  run_cmd "freeview -v OAS2_0002_MR2.long.OAS2_0002/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_MR2.long.OAS2_0002/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR2.long.OAS2_0002/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR2.long.OAS2_0002/surf/rh.white:edgecolor=blue"

  run_cmd "cp OAS2_0002_MR1.long.OAS2_0002/mri/brain.finalsurfs.mgz OAS2_0002_MR1.long.OAS2_0002/mri/brain.finalsurfs.manedit.mgz"


  run_cmd "freeview -v OAS2_0002_MR1.long.OAS2_0002/mri/brain.finalsurfs.manedit.mgz \
            OAS2_0002_MR1.long.OAS2_0002/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR1.long.OAS2_0002/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR1.long.OAS2_0002/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR1.long.OAS2_0002/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR1.long.OAS2_0002/surf/rh.white:edgecolor=blue"

  # run_cmd "mv OAS2_0002_MR2.long.OAS2_0002 OAS2_0002_MR2.long.OAS2_0002_orig"

  # run_cmd "recon-all -long OAS2_0002_MR2 OAS2_0002 -all"

  # run_cmd "mv OAS2_0002_MR1.long.OAS2_0002 OAS2_0002_MR1.long.OAS2_0002_orig"

  # run_cmd "recon-all -long OAS2_0002_MR1 OAS2_0002 -all"

  run_cmd "freeview -v OAS2_0002_MR2.long.OAS2_0002_fixed/mri/brain.finalsurfs.mgz \
            OAS2_0002_MR2.long.OAS2_0002_fixed/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR2.long.OAS2_0002_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR2.long.OAS2_0002_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR2.long.OAS2_0002_fixed/surf/rh.white:edgecolor=blue"

  run_cmd "freeview -v OAS2_0002_MR1.long.OAS2_0002_fixed/mri/brain.finalsurfs.mgz \
            OAS2_0002_MR1.long.OAS2_0002_fixed/mri/aseg.mgz:colormap=lut:opacity=0.25 \
         -f OAS2_0002_MR1.long.OAS2_0002_fixed/surf/lh.pial:edgecolor=red \
            OAS2_0002_MR1.long.OAS2_0002_fixed/surf/rh.pial:edgecolor=red \
            OAS2_0002_MR1.long.OAS2_0002_fixed/surf/lh.white:edgecolor=blue \
            OAS2_0002_MR1.long.OAS2_0002_fixed/surf/rh.white:edgecolor=blue"

fi

if [ $do_ROI_Analysis_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "freeview -v 004/mri/orig.mgz \
  004/mri/aparc+aseg.mgz:colormap=lut:opacity=0.4 \
  -f 004/surf/lh.white:annot=aparc.annot"

  run_cmd "cat $FREESURFER_HOME/FreeSurferColorLUT.txt"

  run_cmd "mri_label2label \
  --srcsubject fsaverage \
  --srcsubject fsaverage \
  --srclabel fsaverage/label/lh.BA45.label \
  --trgsubject 004 \
  --trglabel 004/label/lh.BA45.label \
  --hemi lh \
  --regmethod surface"

  run_cmd "freeview -v 004/mri/orig.mgz"

  run_cmd "freeview -f 004/surf/lh.inflated"

  run_cmd "cd $SUBJECTS_DIR/004/stats"

  run_cmd "cat aseg.stats"

  run_cmd "cd $SUBJECTS_DIR/004/stats"

  run_cmd "cat lh.aparc.stats"

  export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs/group_analysis_tutorial

  run_cmd "cd $SUBJECTS_DIR"

  run_cmd "asegstats2table --subjects 004 021 040 067 080 092 \
  --segno 11 17 18 \
  --tablefile aseg.vol.table"

  run_cmd "asegstats2table \
  --subjects 004 021 040 067 080 092 \
  --segno 11 17 18 \
  --meas mean \
  --tablefile aseg.mean-intensity.table"

  run_cmd "asegstats2table \
  --subjects 004 021 040 067 080 092 \
  --segno 3007 3021 3022 4022 \
  --stats wmparc.stats \
  --tablefile wmparc.vol.table"

  run_cmd "aparcstats2table --hemi lh \
  --subjects 004 021 040 067 080 092 \
  --tablefile lh.aparc.area.table"

  run_cmd "aparcstats2table --hemi lh \
  --subjects 004 021 040 067 080 092 \
  --meas thickness \
  --parc aparc.a2009s \
  --tablefile lh.aparc.a2009.thickness.table"

fi

if [ $do_Diffusion_Processing_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/diffusion_recons
  export TUTORIAL_DIR=$TUTORIAL_DATA/diffusion_tutorial
  subj=Diff001

  ## DO NOT RUN THIS COMMAND. It would take a long time so it has already been run for you. 
  #run_cmd "dt_recon --i $TUTORIAL_DIR/$subj/orig/*-1.dcm 
  #--s $subj --o $TUTORIAL_DIR/$subj/dtrecon"

  run_cmd "freeview -v $TUTORIAL_DIR/$subj/dtrecon/fa.nii \
  $TUTORIAL_DIR/$subj/dtrecon/adc.nii"

  # This command will probably fail because it uses -reg reg.dat instead of lta
  run_cmd "mri_vol2vol --mov $TUTORIAL_DIR/$subj/dtrecon/lowb.nii \
  --targ $SUBJECTS_DIR/$subj/mri/wmparc.mgz \
  --inv --interp nearest --o $SUBJECTS_DIR/$subj/mri/wmparc2diff.mgz \
  --reg $TUTORIAL_DIR/$subj/dtrecon/register.dat --no-save-reg"

  run_cmd "freeview -v $TUTORIAL_DIR/$subj/dtrecon/fa.nii \
  $SUBJECTS_DIR/$subj/mri/aparc+aseg2diff.mgz:colormap=lut"

  run_cmd "mri_mask $TUTORIAL_DIR/$subj/dtrecon/fa.nii \
  $SUBJECTS_DIR/$subj/mri/wmparc2diff.mgz \
  $TUTORIAL_DIR/$subj/dtrecon/fa-masked.mgz"

  run_cmd "freeview -v $TUTORIAL_DIR/$subj/dtrecon/fa.nii \
  $TUTORIAL_DIR/$subj/dtrecon/fa-masked.mgz"

  ## !!! WARNING: DO NOT RUN THIS COMMAND ON THE COURSE LAPTOPS !!!
  #run_cmd "mri_vol2vol --targ $FREESURFER_HOME/subjects/cvs_avg35/mri/norm.mgz \
  #--m3z $SUBJECTS_DIR/$subj/cvs/combined_tocvs_avg35_elreg_afteraseg-norm.m3z \
  #--noDefM3zPath --reg $TUTORIAL_DIR/$subj/dtrecon/register.dat \
  #--mov $TUTORIAL_DIR/$subj/dtrecon/fa-masked.mgz \
  #--o $TUTORIAL_DIR/$subj/dtrecon/fa-masked.ANAT+CVS-to-avg35.mgz \
  #--interp trilin --no-save-reg

  run_cmd "freeview -v $FREESURFER_HOME/subjects/cvs_avg35/mri/norm.mgz \
  $TUTORIAL_DIR/$subj/dtrecon/fa-masked.ANAT+CVS-to-avg35.mgz"

fi

if [ $do_TRACULA_Tutorial -eq 1 ]; then

  export SUBJECTS_DIR=$TUTORIAL_DATA/diffusion_recons

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/dmrirc.tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dwi_motion.txt"

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri"

  run_cmd "freeview $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dtifit_FA.nii.gz"

  run_cmd "cd $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/mni/"

  run_cmd "freeview -v $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz \
  $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/mni/dtifit_FA.bbr.nii.gz"

  run_cmd "freeview -v 
  $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dtifit_FA.nii.gz \
  $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/rh.ilf_AS_avg33_mni_bbr/path.pd.nii.gz:colormap=jet:isosurface=0,0:color='Red':name=rh.ilf \
  $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/path.pd.nii.gz:colormap=jet:isosurface=0,0:color='Red':name=lh.ilf"

  run_cmd "freeview \
  -tv $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/merged_avg33_mni_bbr.mgz \
  -v $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dmri/dtifit_FA.nii.gz"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.overall.txt"

  run_cmd "tractstats2table \
  --inputs $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.overall.txt \
  --overall --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.elmo.2012.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.elmo.2012.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list"

  run_cmd "tractstats2table --load-pathstats-from-file \
  $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list --overall \
  --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "tractstats2table --load-pathstats-from-file \
  $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.list --overall \
  --only-measures FA_Avg --tablefile $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.FA_Avg.table"
  
  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/lh.ilf.All.table"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"

  run_cmd "trac-all -stat -c $TUTORIAL_DATA/diffusion_tutorial/dmrirc.tutorial"

  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"
 
  run_cmd "cat $TUTORIAL_DATA/diffusion_tutorial/elmo.2012/dpath/lh.ilf_AS_avg33_mni_bbr/pathstats.byvoxel.txt"

  run_cmd "freeview -v $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz \
  -w $TUTORIAL_DATA/diffusion_tutorial/stats/*.path.mean.txt"

fi

if [ $do_MultiModal_Integration -eq 1 ]; then

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs"

  run_cmd "cd $SUBJECTS_DIR/multimodal/fmri/fbirn-101"

  run_cmd "freeview -v template.nii \
  $SUBJECTS_DIR/fbirn-anat-101.v4/mri/orig.mgz:visible=0 -f \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/lh.white:edgecolor=green \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/rh.white:edgecolor=green \
  -viewport coronal"

  run_cmd "tkregister2 --mov template.nii --s fbirn-anat-101.v4 \
  --regheader --reg myregister.dat --surf"

  run_cmd "bbregister --mov template.nii --bold \
  --s fbirn-anat-101.v4 \
  --init-rr --lta register.lta"

  run_cmd "cat register.lta"

  run_cmd "freeview -v $SUBJECTS_DIR/fbirn-anat-101.v4/mri/orig.mgz \
  template.nii:reg=bb.register.lta -f \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/lh.white:edgecolor=green \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/rh.white:edgecolor=green \
  -viewport coronal"

  run_cmd "tkregister2 --mov template.nii --reg register.dat --surf"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs"
 
  run_cmd "cd $SUBJECTS_DIR/multimodal/dti"

  run_cmd "freeview -v \
  $SUBJECTS_DIR/M87102113.v4/mri/orig.mgz lowb.nii:reg=dti-analysis/register.lta \
  -f $SUBJECTS_DIR/M87102113.v4/surf/lh.white:edgecolor=green \
  $SUBJECTS_DIR/M87102113.v4/surf/rh.white:edgecolor=green \
  -viewport coronal"

  run_cmd "tkregister2 --mov lowb.nii --reg dti-analysis/register.lta --surf"

  run_cmd "freeview -v $SUBJECTS_DIR/M87102113.v4/mri/brain.mgz \
  $SUBJECTS_DIR/M87102113.v4/mri/orig.mgz \
  $SUBJECTS_DIR/M87102113.v4/mri/wmparc.mgz:colormap=lut:opacity=0.2 \
  fa.nii:reg=dti-analysis/register.lta:colormap=heat:heatscale=.2,.2,1 \
  -colorscale -viewport coronal"

  run_cmd "mri_vol2vol --mov fa.nii \
  --lta dti-analysis/register.lta \
  --fstarg \
  --o fa.anat.mgh"

  run_cmd "mri_segstats \
  --seg $SUBJECTS_DIR/M87102113.v4/mri/wmparc.mgz \
  --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
  --id 251 --id 3021 --id 3024 --id 3030 --id 12 --id 4 \
  --i fa.anat.mgh --sum fa.stats"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs"

  run_cmd "cd $SUBJECTS_DIR/multimodal/fmri/fbirn-101"

  run_cmd "freeview -v $SUBJECTS_DIR/fbirn-anat-101.v4/mri/orig.mgz \
  template.nii:reg=bb.register.lta -f \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/lh.white:edgecolor=green \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/rh.white:edgecolor=green \
  -viewport coronal"

  run_cmd "tkregister2 --mov template.nii \
  --reg bb.register.lta --surf"

  run_cmd "freeview -v $SUBJECTS_DIR/fbirn-anat-101.v4/mri/orig.mgz \
  $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 \
  sig.nii:reg=bb.register.lta:colormap=heat:heatscale=2,2,4 -colorscale \
  -viewport coronal"

  run_cmd "mri_vol2surf --mov sig.nii \
    --reg bb.register.lta \
    --projfrac 0.5 --interp nearest \
    --hemi lh --o lh.sig.mgh"

  run_cmd "mri_info lh.sig.mgh"

  run_cmd "freeview -f $SUBJECTS_DIR/fbirn-anat-101.v4/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay=lh.sig.mgh:overlay_threshold=2,5 \
-viewport 3d"

  run_cmd "mri_vol2vol --mov ces.nii \
    --reg bb.register.lta \
    --fstarg --interp nearest \
    --o ces.anat.bb.mgh"

  run_cmd "mri_info ces.anat.bb.mgh"

  run_cmd "mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.bb.stats"

  run_cmd "mri_vol2vol --mov sig.nii \
    --reg bb.register.lta \
    --fstarg --interp nearest \
    --o sig.anat.bb.mgh"

  run_cmd "mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.abs-masked.bb.stats \
   --mask sig.anat.bb.mgh --maskthresh 2 --masksign abs"

  run_cmd "mri_segstats \
   --seg $SUBJECTS_DIR/fbirn-anat-101.v4/mri/aparc+aseg.mgz \
   --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt \
   --id 1021 --id 1022 --id 1030  --id 17 \
   --i ces.anat.bb.mgh --sum ces.pos-masked.bb.stats \
   --mask sig.anat.bb.mgh --maskthresh 2 --masksign pos"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/buckner_data/tutorial_subjs"

  run_cmd "cd $SUBJECTS_DIR/multimodal/fmri"

  run_cmd "freeview -v $SUBJECTS_DIR/fbirn-anat-101.v4/mri/orig.mgz \
  fbirn-101/template.nii:reg=fbirn-101/bb.register.lta -f \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/lh.white:edgecolor=green \
  $SUBJECTS_DIR/fbirn-anat-101.v4/surf/rh.white:edgecolor=green \
  -viewport coronal"

  run_cmd "tkregister2 --mov fbirn-101/template.nii \
  --reg fbirn-101/bb.register.lta --surf"

  run_cmd "mris_preproc --target fsaverage --hemi lh \
  --iv  fbirn-101/ces.nii fbirn-101/bb.register.lta \
  --iv  fbirn-103/ces.nii fbirn-103/bb.register.lta \
  --iv  fbirn-104/ces.nii fbirn-104/bb.register.lta \
  --iv  fbirn-105/ces.nii fbirn-105/bb.register.lta \
  --iv  fbirn-106/ces.nii fbirn-106/bb.register.lta \
  --projfrac 0.5 \
  --out lh.ces.mgh"

  run_cmd "mri_info lh.ces.mgh"

  run_cmd "mri_surf2surf --hemi lh --s fsaverage --fwhm 5 --cortex\
  --sval lh.ces.mgh --tval lh.ces.sm05.mgh"

  run_cmd "mri_glmfit --y lh.ces.sm05.mgh --surf fsaverage lh \
  --osgm --glmdir lh.ces.sm05.osgm --cortex"

  run_cmd "freeview -f \
  $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=aparc.annot:annot_outline=1:overlay=lh.ces.sm05.osgm/osgm/sig.mgh:overlay_threshold=2,5 \
  -viewport 3d"

  run_cmd "asegstats2table \
  --meas volume \
  --tablefile ces.pos-masked.vol.stats \
  --i fbirn-101/ces.pos-masked.bb.stats \
  fbirn-103/ces.pos-masked.bb.stats \
  fbirn-104/ces.pos-masked.bb.stats \
  fbirn-105/ces.pos-masked.bb.stats \
  fbirn-106/ces.pos-masked.bb.stats"

  run_cmd "asegstats2table \
  --meas mean \
  --tablefile ces.abs-masked.mean.stats \
  --i fbirn-101/ces.abs-masked.bb.stats \
  fbirn-103/ces.abs-masked.bb.stats \
  fbirn-104/ces.abs-masked.bb.stats \
  fbirn-105/ces.abs-masked.bb.stats \
  fbirn-106/ces.abs-masked.bb.stats"

  run_cmd "asegstats2table \
  --meas mean \
  --tablefile ces.pos-masked.mean.stats \
  --i fbirn-101/ces.pos-masked.bb.stats \
  fbirn-103/ces.pos-masked.bb.stats \
  fbirn-104/ces.pos-masked.bb.stats \
  fbirn-105/ces.pos-masked.bb.stats \
  fbirn-106/ces.pos-masked.bb.stats"

fi

if [ $do_FSFAST_Tutorial -eq 1 ]; then

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"
 
  run_cmd "ls"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/sess01.noproc"

  run_cmd "ls"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/fsfast-tutorial.subjects"

  run_cmd "cat subjectname"

  run_cmd "ls $SUBJECTS_DIR"

  run_cmd "ls rest"

  run_cmd "ls bold"

  run_cmd "cd bold/001"

  run_cmd "ls"

  run_cmd "mri_info --dim f.nii.gz"

  run_cmd "mri_info --res f.nii.gz"

  run_cmd "freeview -v f.nii.gz -timecourse"

  run_cmd "cat workmem.par"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"

  run_cmd "cat sessidlist"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/fsfast-tutorial.subjects"
  
  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"

  run_cmd "preproc-sess -s sess01 -fsd bold \
  -stc up -surface fsaverage lhrh -mni305 -fwhm 5 -per-run"

  run_cmd "ls $TUTORIAL_DATA/fsfast-functional/sess01/bold/001"

  #run_cmd "plot-twf-sess -s sess01 -fsd bold -mc"

  run_cmd "tkregister-sess -s sess01 -s sess02 \
  -s sess03 -fsd bold -per-run -bbr-sum"

  run_cmd "tkregister-sess -s sess02 -fsd bold -per-run"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/sess01/bold/001"

  run_cmd "ls -ltr"

  run_cmd "freeview -v template.nii.gz \
  masks/brain.nii.gz:colormap=heat -viewport coronal"

  run_cmd "freeview -v template.nii.gz \
  masks/brain.e3.nii.gz:colormap=heat -viewport coronal"

  run_cmd "mri_info --dim fmcpr.up.sm5.mni305.2mm.nii.gz"

  run_cmd "mri_info --res fmcpr.up.sm5.mni305.2mm.nii.gz"

  run_cmd "mri_info --dim fmcpr.up.sm5.fsaverage.lh.nii.gz"

  run_cmd "mri_info --res fmcpr.up.sm5.fsaverage.lh.nii.gz"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/fsfast-tutorial.subjects"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"

  run_cmd "mkanalysis-sess \
  -fsd bold -stc up  -surface fsaverage lh -fwhm 5  \
  -event-related  -paradigm workmem.par -nconditions 5 \
  -spmhrf 0 -TR 2 -refeventdur 16 -nskip 4 -polyfit 2 \
  -analysis my-workmem.sm05.lh -force -per-run"

  run_cmd "ls my-workmem.sm05.lh | grep analysis.info"

  run_cmd "mkcontrast-sess -analysis my-workmem.sm05.lh -contrast encode-v-base -a 1"
  run_cmd "ls my-workmem.sm05.lh"

  run_cmd "mkcontrast-sess -analysis my-workmem.sm05.lh -contrast emot.dist-v-neut.dist -a 2 -c 3"

  run_cmd "mkcontrast-sess -analysis my-workmem.sm05.lh -contrast distractor.avg-v-base -a 2 -a 3"

  run_cmd "mkcontrast-sess -analysis my-workmem.sm05.lh -contrast emot.dist-v-base -a 2"

  run_cmd "mkcontrast-sess -analysis my-workmem.sm05.lh -contrast probe.avg-v-base -a 4 -a 5"

  run_cmd "mkanalysis-sess \
  -fsd bold -stc up  -surface fsaverage rh -fwhm 5  \
  -event-related  -paradigm workmem.par -nconditions 5 \
  -spmhrf 0 -TR 2 -refeventdur 16 -nskip 4 -polyfit 2 \
  -analysis my-workmem.sm05.rh -force -per-run"

  run_cmd "mkanalysis-sess \
  -fsd bold -stc up  -mni305 2 -fwhm 5  \
  -event-related  -paradigm workmem.par -nconditions 5 \
  -spmhrf 0 -TR 2 -refeventdur 16 -nskip 4 -polyfit 2 \
  -analysis my-workmem.sm05.mni305 -force -per-run"

  run_cmd "selxavg3-sess -s sess01 -analysis workmem.sm05.lh"

  run_cmd "ls $TUTORIAL_DATA/fsfast-functional/sess01/bold"

  run_cmd "tksurfer-sess -s sess01 \
  -analysis workmem.sm05.lh \
  -c encode-v-base \
  -c emot.dist-v-base \
  -c probe.avg-v-base \
  -c emot.dist-v-neut.dist \
  -freeview"

  run_cmd "tksurfer-sess -s sess01 \
  -analysis workmem.sm05.rh \
  -c encode-v-base \
  -c emot.dist-v-base \
  -c emot.dist-v-neut.dist \
  -c probe.avg-v-base \
  -freeview"

  run_cmd "tkmedit-sess -s sess01 -analysis workmem.sm05.mni305 \
  -c encode-v-base \
  -c emot.dist-v-base \
  -c emot.dist-v-neut.dist \
  -c probe.avg-v-base \
  -freeview"

  run_cmd "cd sess01/bold"
 
  run_cmd "ls | grep workmem.sm05.lh"

  run_cmd "cd workmem.sm05.lh"

  run_cmd "ls"

  run_cmd "cd encode-v-base"

  run_cmd "ls"

  run_cmd "export SUBJECTS_DIR=$TUTORIAL_DATA/fsfast-tutorial.subjects"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional"

  run_cmd "isxconcat-sess -sf sessidlist -analysis workmem.sm05.lh -contrast encode-v-base -o my-group"

  run_cmd "isxconcat-sess -sf sessidlist -analysis workmem.sm05.rh -contrast encode-v-base -o my-group"

  run_cmd "isxconcat-sess -sf sessidlist -analysis workmem.sm05.mni305 -contrast encode-v-base -o my-group"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/my-group"

  run_cmd "ls"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/my-group/workmem.sm05.lh"

  run_cmd "ls"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/group/workmem.sm05.lh/encode-v-base"

  run_cmd "ls"

  run_cmd "mri_glmfit --y ces.nii.gz \
  --wls cesvar.nii.gz \
  --osgm \
  --surface fsaverage lh \
  --glmdir my-glm.wls \
  --nii.gz"

  run_cmd "tksurferfv fsaverage lh inflated -aparc -overlay my-glm.wls/osgm/sig.nii.gz -fminmax 2 3"

  run_cmd "mri_glmfit-sim --glmdir my-glm.wls --cache 3 pos --cwpvalthresh .0166"

  run_cmd "cat my-glm.wls/osgm/cache.th30.pos.sig.cluster.summary"

  run_cmd "tksurferfv fsaverage lh inflated \
  -overlay my-glm.wls/osgm/cache.th30.pos.sig.cluster.nii.gz \
  -annot ./my-glm.wls/osgm/cache.th30.pos.sig.ocn.annot -fminmax 1.3 3"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/group/workmem.sm05.rh/encode-v-base"

  run_cmd "mri_glmfit --y ces.nii.gz --wls cesvar.nii.gz --osgm \
  --surface fsaverage rh --glmdir my-glm.wls --nii.gz"

  run_cmd "mri_glmfit-sim --glmdir my-glm.wls --cache 3 pos --cwpvalthresh .0166"

  run_cmd "cat my-glm.wls/osgm/cache.th30.pos.sig.cluster.summary"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/group/workmem.sm05.mni305/encode-v-base" 

  run_cmd "ls"
  
  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/group/workmem.sm05.mni305/encode-v-base"

  run_cmd "mri_glmfit --y ces.nii.gz --wls cesvar.nii.gz --osgm  \
  --glmdir my-glm.wls --nii.gz"

  run_cmd "freeview -v $SUBJECTS_DIR/fsaverage/mri/orig.mgz $SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 my-glm.wls/osgm/sig.nii.gz:colormap=heat:heatscale=2,2,3"

  run_cmd "mri_glmfit-sim --glmdir my-glm.wls --grf 3 pos --cwpvalthresh .0166"

  run_cmd "cat my-glm.wls/osgm/grf.th3.pos.sig.cluster.summary"

  run_cmd "freeview $SUBJECTS_DIR/fsaverage/mri/orig.mgz \
my-glm.wls/osgm/grf.th3.pos.sig.ocn.anat.nii.gz:colormap=lut:lut=./my-glm.wls/osgm/grf.th3.pos.sig.ocn.lut:opacity=0 \
my-glm.wls/osgm/grf.th3.pos.sig.cluster.nii.gz:colormap=heat:heatscale=1.3,1.3,5"

  run_cmd "cd $TUTORIAL_DATA/fsfast-functional/group"

  run_cmd "cat workmem.sm05.lh/encode-v-base/my-glm.wls/osgm/cache.th30.pos.sig.cluster.summary"

  run_cmd "cat workmem.sm05.rh/encode-v-base/my-glm.wls/osgm/cache.th30.pos.sig.cluster.summary"

  run_cmd "cat workmem.sm05.mni305/encode-v-base/my-glm.wls/osgm/grf.th3.pos.sig.cluster.summary"

  run_cmd "vlrmerge --o encode.merged.nii.gz \
  --lh workmem.sm05.lh/encode-v-base/my-glm.wls/osgm/cache.th30.pos.sig.cluster.nii.gz \
  --rh workmem.sm05.rh/encode-v-base/my-glm.wls/osgm/cache.th30.pos.sig.cluster.nii.gz \
  --vol workmem.sm05.mni305/encode-v-base/my-glm.wls/osgm/grf.th3.pos.sig.cluster.nii.gz \
  --scm workmem.sm05.mni305/subcort.mask.nii.gz"

  run_cmd "freeview -v $SUBJECTS_DIR/fsaverage/mri/orig.mgz \
  $SUBJECTS_DIR/fsaverage/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 \
  encode.merged.nii.gz:colormap=heat:heatscale=1.3,1.3,5 \
  -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:edgecolor=red \
  $SUBJECTS_DIR/fsaverage/surf/rh.pial:edgecolor=red \
  $SUBJECTS_DIR/fsaverage/surf/lh.white:edgecolor=yellow \
  $SUBJECTS_DIR/fsaverage/surf/rh.white:edgecolor=yellow"

fi
