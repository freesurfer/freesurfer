##
## Makefile.am
##

ACLOCAL_AMFLAGS=-I m4

foodir=$(prefix)
foo_DATA=LICENSE NOTICE README AUTHORS INSTALL COPYING

OCTSUBDIRS=\
	oct_train \
	oct_register_mosaic 

DMRISUBDIRS=\
	dmri_poistats \
	dmri_tensoreig	

GLSUBDIRS=\
	glut \
	nmovie \
	mris_show \
	stim_polar \
	tkmedit \
	tkregister2 \
	tksurfer \
	mri_probedicom \
	orient_mri \
	qdec \
	mris2rgb

MRISUBDIRS=\
	check_siemens_dir \
	repair_siemens_file \
	mkxsubjreg \
	rbftest \
	lta_convert \
	mri_interpolate \
	mri_elastic_energy \
	mri_jacobian \
	mri_update_gca \
	mri_coreg \
	mri_compile_edits \
	mri_compute_volume_fractions \
	mri_label_accuracy \
	mri_add_xform_to_header \
	mri_bias \
	mri_annotation2label \
	mri_aparc2aseg \
	mri_aparc2wmseg \
	mri_auto_fill \
	mri_average \
	mri_dct_align \
	mri_binarize \
	mri_build_priors \
	mri_ca_label \
	mri_correct_segmentations \
	mri_cal_renormalize_gca \
	mri_divide_segmentation \
	mri_gtmpvc \
	mri_rf_label \
	mri_rf_train \
	mri_rf_long_label \
	mri_rf_long_train \
	mri_aseg_edit_train \
	mri_aseg_edit_reclassify \
	mri_ca_normalize \
	mri_ca_tissue_parms \
	mri_ca_register \
	mri_ca_train \
	mri_cvs_register \
	mri_gcab_train \
	mri_gtmseg \
	mri_gcut \
	mri_cc mri_cht2p \
	mri_cnr \
	mri_compute_overlap \
	mri_concat \
	mri_concatenate_gcam \
	mri_convert \
	mri_mosaic \
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
	mri_fcili \
	mri_fdr \
	mri_fieldsign \
	mri_fill \
	mri_fit_bias \
	mri_fwhm \
	mri_gca_ambiguous \
	mri_head \
	histo_segment \
	histo_synthesize \
	histo_fix_topology \
	mri_histo_eq \
	mri_compute_structure_transforms \
	mri_ibmc \
	mri_iic \
	mri_info \
	mri_and \
	mri_hausdorff_dist \
	mri_histo_normalize \
	mri_fuse_segmentations \
	mri_fuse_intensity_images \
	mri_simulate_atrophy \
	mri_joint_density \
	mri_label_histo \
	mri_label_vals \
	mri_label_volume \
	mri_label2label \
	mri_label2vol \
	mri_linear_register \
	mri_make_bem_surfaces \
	mri_make_labels \
	mri_make_register \
	mri_make_template \
	mri_map_cpdat \
	mri_mark_temporal_lobe \
	mri_mask mri_ms_fitparms \
	mri_mcsim \
	mri_mc \
	mri_matrix_multiply \
	mri_modify \
	mri_multispectral_segment \
	mri_multiscale_segment \
	mri_nlfilter \
	mri_normalize \
	mri_otl \
	mri_partial_ribbon \
	mri_paint \
	mri_parselabel \
	mri_parse_sdcmdir \
	mri_polv \
	mri_probe_ima \
	mri_reduce \
	mri_register \
	mri_relabel_hypointensities \
	mri_relabel_nonwm_hypos \
	mri_remove_neck \
	mri_ribbon \
	mri_rigid_register \
	mri_robust_register \
	mri_seg_diff \
	mri_segment \
	mri_segment_tumor \
	mri_seghead \
	mri_surfacemask \
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
	mri_train_autoencoder \
	mri_classify \
	mri_twoclass \
	mri_volcluster \
	mri_voldiff \
	mri_volsynth \
	mri_vol2roi \
	mri_vol2vol \
	mri_segreg \
	mri_segstats \
	mri_stats2seg \
	mri_vol2surf \
	mri_warp_convert \
	mri_watershed \
	mri_wbc \
	mri_wmfilter \
	mri_xcorr \
	mri_xvolavg \
	mri_z2p \
	mri_threshold \
	mri_edit_wm_with_aseg \
	mri_bc_sc_bias_correct \
	mri_morphology \
	mri_topologycorrection \
	mri_segment_wm_damage \
	mri_make_density_map \
	mri_path2label \
	mri_fslmat_to_lta \
	mri_ms_EM \
	mri_ms_EM_with_atlas \
	mri_extract_largest_CC \
	mri_compute_seg_overlap \
	mri_concatenate_lta \
	mri_normalize_tp2 \
	mri_ms_LDA \
	mri_transform_to_COR \
	mri_log_likelihood \
	mri_compute_change_map \
	mri_apply_inu_correction \
	talairach_afd \
	talairach_avi \
	template

MRISSUBDIRS=\
	label_area \
	label_border \
	label2patch \
	label2flat \
	spherical_stats \
	tridec \
	mris_copy_header \
	mris_segmentation_stats \
	mris_simulate_atrophy \
	mris_register_to_volume \
	mris_init_global_tractography \
	mris_interpolate_warp \
	mris_parcellate_connectivity \
	mris_map_cuts \
	mris_hausdorff_dist \
	mris_merge_parcellations \
	mris_distance_transform \
	mris_svm_classify \
	mris_svm_train \
	mris_calc \
	mris_shrinkwrap \
	mris_density \
	mris_rf_label \
	mris_rf_train \
	mris_fbirn_annot \
	mris_fill \
	mris_add_template \
	mris_anatomical_stats \
	mris_transmantle_dysplasia_paths \
	mris_annot_to_segmentation \
	mris_apply_reg \
	mris_aseg_distance \
	mris_average_curvature \
	mris_ca_label \
	mris_divide_parcellation \
	mris_ca_train \
	mris_classify_thickness \
	mris_compute_acorr \
	mris_compute_lgi \
	mris_compute_overlap \
	mris_compute_parc_overlap \
	mris_convert \
	mris_decimate \
	mris_decimate_gui \
	mris_distance_map \
	mris_curvature \
	mris_curvature_stats \
	mris_entropy \
	mris_errors \
	mris_euler_number \
	mris_expand \
	mris_find_flat_regions \
	mris_fix_topology \
	mris_topo_fixer \
	mris_flatten \
	mris_fwhm \
	mris_inflate \
	mris_info \
	mris_jacobian \
	mris_label_area \
	mris_label_calc \
	mris_label_mode \
	mris_label2annot \
	mris_make_average_surface \
	mris_make_surfaces \
	mris_make_map_surfaces \
	mris_longitudinal_surfaces \
	mris_make_template \
	mris_morph_stats \
	mris_deform \
	mris_multiscale_stats \
	mris_niters2fwhm \
	mris_pmake \
	mris_remove_variance \
	mris_rescale \
	mris_reverse \
	mris_rotate \
	mris_ms_refine \
	mris_diff \
	mris_distance_to_label \
	mris_register_label_map \
	mris_register \
	mris_left_right_register \
	mris_congeal \
	mris_resample \
	mris_sample_label \
	mris_sample_parc \
	mris_segment_vals \
	mris_segment \
	mris_seg2annot \
	mris_smooth \
	mris_make_face_parcellation \
	mris_sphere \
	mris_spherical_average \
	mris_average_parcellation \
	mris_surface_to_vol_distances \
	mrisp_paint \
	mris_talairach \
	mris_surface_change \
	mris_thickness \
	mris_BA_segment \
	mris_compute_optimal_kernel \
	mris_thickness_comparison \
	mris_transform \
	mris_translate_annotation \
	mris_twoclass \
	mris_w_to_curv \
	mris_volmask \
	mris_volume \
	mris_warp \
	mris_thickness_diff \
	mris_ms_surface_CNR \
	mris_mesh_subdivide \
	mris_surface_stats \
	mris_compute_volume_fractions

FSGDF=\
	fsgdf \
	mri_gdfglm \
	mri_glmfit \
	mris_glm

FSLUTILS=fslutils

STATSDIRS=\
	optseq2 \
	stat_normalize

MISCDIRS=\
	histo_compute_joint_density \
	histo_register_block \
	regdat2xfm \
	test_window_env \
	dngtester \
	test_makevol \
	images \
	dummy \
	dummy_qt \
	fem_elastic \
	hiam_make_surfaces\
	hiam_make_template\
	hiam_register \
	histo_register \
	spline3 \
	swi_processing \
	nmovie_qt \
	HippoSF \
	BrainstemSS

# ordering is important here for proper 'make check' and 'make distcheck' function
COREDIRS=\
	include \
	unix \
	dicom \
	hipsstubs \
	expat \
	jpeg \
	tiff \
	xml2 \
	rgb \
	log \
	gpu_utils \
	opencl_algorithms \
	utils \
	utilscpp \
	vtkfsio \
	vtkutils \
	itkio \
	itkutils \
	tetgen \
	lineprof \
	$(FSGDF) \
	connectgraph \
	qdecproject \
	svm-lib \
	scripts \
	distribution \
	$(FSLUTILS) \
	$(STATSDIRS) \
	$(GLSUBDIRS)

# dependent on COREDIRS
OTHERDIRS=\
	$(MRISUBDIRS) \
	$(MRISSUBDIRS) \
	$(MISCDIRS) \
	$(DMRISUBDIRS) \
	$(OCTSUBDIRS) \
	fsfast \
	matlab \
	qdec_glmfit \
	cudatest

# these have lots of files, each will maximize the multi-threaded make option (-j)
BIGDIRS=\
	freeview \
	GEMS \
	trc

# SUBDIRS is the automake definition of the 'all' for make
SUBDIRS=\
	$(COREDIRS) \
	$(OTHERDIRS) \
	$(BIGDIRS)

# this target is an ad-hoc way to run the build a lot faster
# by launching makes in parallel, after the dependencies are
# built.  also make is run multi-threaded, which speeds the
# build within a directory.
fastbuild:
	@list='$(COREDIRS)'; \
	for d in $$list; do \
		cd "$$d"; make -j 16; cd - ; \
	done
	@list='$(OTHERDIRS)'; \
	for d in $$list; do \
		cd "$$d"; (make -j 16 &) ; sleep 0.4; cd - ; \
	done
	@list='$(BIGDIRS)'; \
	for d in $$list; do \
		cd "$$d"; (make -j 16); cd - ; \
	done
	sleep 10
	make

# build stamp and symlinks to necessary 3rd-party packages.
# in the case of FSL, just copy the 'flirt' utility to our bin
install-data-hook:
	echo "$(PACKAGE)-$(target_triplet)-$(VERSION)-`date +%Y%m%d`" \
		> $(DESTDIR)$(prefix)/build-stamp.txt ;
	cp -f $(DESTDIR)$(prefix)/build-stamp.txt \
		$(DESTDIR)$(prefix)/VERSION
	if test ! "x$(MNI_DIR)" = "x" ; then \
		if test ! "$(MNI_DIR)" = "$(DESTDIR)$(prefix)/mni" ; then \
			if test -e $(MNI_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/mni ; \
				ln -s $(MNI_DIR) $(DESTDIR)$(prefix)/mni ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(FSLDIR)" = "x" ; then \
		if test -e $(FSLDIR) ; then \
			cp -f $(FSLDIR)/bin/flirt $(DESTDIR)$(bindir)/flirt.fsl ; \
			cp -f $(FSLDIR)/bin/fslmaths $(DESTDIR)$(bindir)/fslmaths.fsl ; \
			cp -f $(FSLDIR)/bin/slicetimer $(DESTDIR)$(bindir)/slicetimer.fsl ; \
			cp -f $(FSLDIR)/bin/bet $(DESTDIR)$(bindir)/bet.fsl ; \
			cp -f $(FSLDIR)/bin/fslswapdim $(DESTDIR)$(bindir)/fslswapdim.fsl ; \
			cp -f $(FSLDIR)/bin/fslorient $(DESTDIR)$(bindir)/fslorient.fsl ; \
		fi; \
	fi;
	if test ! "x$(DCMTKDIR)" = "x" ; then \
		if test -e $(DCMTKDIR) ; then \
			cp -f $(DCMTKDIR)/bin/dcmdjpeg $(DESTDIR)$(bindir)/dcmdjpeg.fs ; \
			cp -f $(DCMTKDIR)/bin/dcmdrle $(DESTDIR)$(bindir)/dcmdrle.fs ; \
		fi; \
	fi;
	if test ! "x$(GDCMDIR)" = "x" ; then \
		if test -e $(GDCMDIR) ; then \
			cp -f $(GDCMDIR)/bin/gdcmconv $(DESTDIR)$(bindir)/gdcmconv.fs ; \
		fi; \
	fi;
	if test ! "x$(AFNIDIR)" = "x" ; then \
		if test -e $(AFNIDIR)/bin/3dvolreg ; then \
			cp -f $(AFNIDIR)/bin/3dvolreg $(DESTDIR)$(bindir)/3dvolreg.afni ; \
		fi; \
		if test -e $(AFNIDIR)/3dvolreg ; then \
			cp -f $(AFNIDIR)/3dvolreg $(DESTDIR)$(bindir)/3dvolreg.afni ; \
		fi; \
	fi;
	if test ! "x$(TCL_DIR)" = "x" ; then \
		if test ! "$(TCL_DIR)" = "$(DESTDIR)$(prefix)/lib/tcltktixblt" ; then \
			if test -e $(TCL_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/tcltktixblt ; \
				ln -s $(TCL_DIR) $(DESTDIR)$(prefix)/lib/tcltktixblt ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(VTK_DIR)" = "x" ; then \
		if test ! "$(VTK_DIR)" = "$(DESTDIR)$(prefix)/lib/vtk" ; then \
			if test -e $(VTK_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/vtk ; \
				ln -s $(VTK_DIR) $(DESTDIR)$(prefix)/lib/vtk ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(QT_DIR)" = "x" ; then \
		if test ! "$(QT_DIR)" = "$(DESTDIR)$(prefix)/lib/qt" ; then \
			if test -e $(QT_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/qt ; \
				ln -s $(QT_DIR) $(DESTDIR)$(prefix)/lib/qt ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(KWWIDGETS_DIR)" = "x" ; then \
		if test ! "$(KWWIDGETS_DIR)" = "$(DESTDIR)$(prefix)/lib/KWWidgets" ; then \
			if test -e $(KWWIDGETS_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/KWWidgets ; \
				ln -s $(KWWIDGETS_DIR) $(DESTDIR)$(prefix)/lib/KWWidgets ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(PETSC_DIR)" = "x" ; then \
		if test ! "$(PETSC_DIR)" = "$(DESTDIR)$(prefix)/lib/petsc" ; then \
			if test -e $(PETSC_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/petsc ; \
				ln -s $(PETSC_DIR) $(DESTDIR)$(prefix)/lib/petsc ; \
			fi; \
		fi; \
	fi;
	if test ! "x$(CUDA_DIR)" = "x" ; then \
		if test ! "$(CUDA_DIR)" = "$(DESTDIR)$(prefix)/lib/cuda" ; then \
			if test -e $(CUDA_DIR) ; then \
				rm -Rf $(DESTDIR)$(prefix)/lib/cuda ; \
				ln -s $(CUDA_DIR) $(DESTDIR)$(prefix)/lib/cuda ; \
			fi; \
		fi; \
	fi;

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/build-stamp.txt
	rm -f $(DESTDIR)$(prefix)/VERSION
	rm -f $(DESTDIR)$(bindir)/3dvolreg.afni
	rm -f $(DESTDIR)$(bindir)/flirt.fsl
	rm -f $(DESTDIR)$(bindir)/slicetimer.fsl
	rm -f $(DESTDIR)$(bindir)/bet.fsl

EXTRA_DIST=setup_configure Makefile.extra x86cpucaps \
	doxy-file-header Doxyfile $(foo_DATA)

DISTCLEANFILES=cpuflags.mak

# we have files in our dist that dont and shouldnt get cleaned, like testdata.
# this target override allows make distcheck to pass
distcleancheck:
	echo "Skipping target distcleancheck"

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=CVS core
include $(top_srcdir)/Makefile.extra
