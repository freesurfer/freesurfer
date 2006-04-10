##
## Makefile.am
##

foodir=$(prefix)
foo_DATA=README COPYING

GLSUBDIRS=mris2rgb \
	nmovie \
	mris_show \
	tkmedit \
	tkregister \
	tkregister2 \
	tksurfer \
	mri_probedicom

MRISUBDIRS=check_siemens_dir \
	repair_siemens_file \
	mkxsubjreg \
	rbftest \
	rawtocor \
	mri_add_xform_to_header \
	mri_annotation2label \
	mri_aparc2aseg \
	mri_aparc2wmseg \
	mri_auto_fill \
	mri_average \
	mri_build_priors \
	mri_ca_label \
	mri_ca_normalize \
	mri_ca_tissue_parms \
	mri_ca_register \
	mri_ca_train \
	mri_cc mri_cht2p \
	mri_cnr \
	mri_compute_overlap \
	mri_concat \
	mri_convert \
	mri_convert_mdh \
	mri_copy_values \
	mri_cor2label \
	mri_deface \
	mri_diff \
	mri_hires_register \
	mri_edit_segmentation \
	mri_edit_segmentation_with_surfaces \
	mri_em_register \
	mri_estimate_tissue_parms \
	mri_evaluate_morph \
	mri_extract \
	mri_extract_conditions \
	mri_extract_label \
	mri_fill \
	mri_fwhm \
	mri_gca_ambiguous \
	mri_glmfit \
	mri_head \
	mri_histo_eq \
	mri_info \
	mri_joint_density \
	mri_label_histo \
	mri_label_vals \
	mri_label_volume \
	mri_label2label \
	mri_label2vol \
	mri_linear_register \
	mri_make_bem_surfaces \
	mri_make_register \
	mri_make_template \
	mri_mark_temporal_lobe \
	mri_mask mri_ms_fitparms \
	mri_mcsim \
	mri_mc \
	mri_matrix_multiply \
	mri_modify \
	mri_nlfilter \
	mri_normalize \
	mri_partial_ribbon \
	mri_paint \
	mri_parselabel \
	mri_parse_sdcmdir \
	mri_polv \
	mri_probe_ima \
	mri_reduce \
	mri_register \
	mri_relabel_hypointensities \
	mri_remove_neck \
	mri_ribbon \
	mri_rigid_register \
	mri_segment \
	mri_seghead \
	mri_surfacemask \
	mri_strip_skull \
	mri_strip_subject_info \
	mri_strip_nonwhite \
	mri_surf2surf \
	mri_surf2vol \
	mri_surfcluster \
	mri_synthesize \
	mri_tessellate \
	mri_distance_transform \
	mri_transform \
	mri_train \
	mri_twoclass \
	mri_volcluster \
	mri_voldiff \
	mri_volsynth \
	mri_vol2roi \
	mri_vol2vol \
	mri_segstats \
        mri_stats2seg \
	mri_vol2surf \
	mri_watershed \
	mri_wmfilter \
	mri_xvolavg \
	mri_threshold \
	mri_edit_wm_with_aseg \
	mri_bc_sc_bias_correct \
	mri_morphology \
	mri_topologycorrection \
	mri_segment_wm_damage \
	mri_make_density_map \
	mri_path2label

MRISSUBDIRS=stim \
	stim_polar \
	tridec \
	mris_add_template \
	mris_anatomical_stats \
	mris_annot_to_segmentation \
	mris_average_curvature \
	mris_ca_label \
	mris_ca_train \
	mris_classify_thickness \
	mris_compute_acorr \
	mris_compute_overlap \
	mris_convert \
	mris_curvature \
	mris_curvature_stats \
	mris_entropy \
	mris_errors \
	mris_euler_number \
	mris_expand \
	mris_find_flat_regions \
	mris_fix_topology \
	mris_flatten \
	mris_fwhm \
	mris_inflate \
	mris_info \
	mris_jacobian \
	mris_label_area \
	mris_label_mode \
	label_area \
	mris_make_average_surface \
	mris_make_surfaces \
	mris_longitudinal_surfaces \
	mri_multispectral_segment \
	mris_make_template \
	mris_morph_stats \
	mris_multiscale_stats \
	mris_niters2fwhm \
	mris_remove_variance \
	mris_rescale \
	mris_reverse \
	mris_rotate \
	mris_ms_refine \
	mris_diff \
	mris_distance_to_label \
	mris_register \
	mris_sample_label \
	mris_sample_parc \
	mris_segment_vals \
	mris_smooth mris_sphere \
	mris_spherical_average \
	mris_surface_to_vol_distances \
	mrisp_paint \
	mris_talairach \
	mris_thickness \
	mris_thickness_comparison \
	mris_transform \
	mris_translate_annotation \
	mris_twoclass \
	mris_w_to_curv \
	spherical_stats \
	mris_volume

FSGDF=	mri_gdfglm \
	mris_glm

FSLUTILS=fslutils

STATSDIRS=optseq2 \
	spm2bfloat \
	statcombine \
	stat_normalize

MISCDIRS=histo_compute_joint_density \
	histo_register_block \
	regdat2xfm \
	brik2bfloat \
	minc_to_bshort \
	test_window_env \
	dngtester \
	test_makevol \
	scripts \
	images \
	distribution

SUBDIRS=unix \
	dicom \
	hipsstubs \
	rgb \
	utils \
	utilscpp \
	log \
	fsgdf \
	$(MRISUBDIRS) \
	$(MRISSUBDIRS) \
	$(GLSUBDIRS) \
	$(FSGDF) \
	$(FSLUTILS) \
	$(STATSDIRS) \
	$(MISCDIRS) \
	dummy \
	scuba \
	fsfast \
	matlab \
	hiam_make_surfaces \
	include \
	qdec \
	mri_fslmat_to_lta \
	plotter \
	mri_ms_EM \
	mri_ms_EM_with_atlas \
	mris_thickness_diff \
	mri_extract_largest_CC \
	mri_compute_seg_overlap \
	mris_ms_surface_CNR \
	mri_concatenate_lta \
	mris_surface_stats \
	mri_normalize_tp2 \
	mri_ms_LDA \
	mri_transform_to_COR \
	mri_apply_inu_correction

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
