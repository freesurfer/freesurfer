
set (attic_binaries "${CMAKE_INSTALL_PREFIX} check_siemens_dir connectgraph dmri_tensoreig histo_fix_topology histo_segment label_area label_border mri_apply_inu_correction mri_aseg_edit_reclassify mri_aseg_edit_train mri_auto_fill mri_bc_sc_bias_correct mri_build_priors mri_cht2p mri_classify mri_compute_structure_transforms mri_concatenate_gcam mri_convert_mdh mri_divide_segmentation mri_elastic_energy mri_estimate_tissue_parms mri_extract_conditions mri_fcili mri_gcab_train mri_gdfglm mri_histo_normalize mri_ibmc mri_interpolate mri_label_accuracy mri_make_bem_surfaces mri_make_density_map mri_make_labels mri_make_register mri_make_template mri_mosaic mri_ms_EM mri_ms_EM_with_atlas mri_ms_LDA mri_multiscale_segment mri_multispectral_segment mri_parselabel mri_partial_ribbon mri_register mri_rf_long_label mris_add_template mris_aseg_distance mris_average_parcellation mris_classify_thickness mris_compute_optimal_kernel mris_density mri_segment_tumor mri_segment_wm_damage mris_fbirn_annot mris_glm mri_simulate_atrophy mris_longitudinal_surfaces mris_make_map_surfaces mris_ms_surface_CNR mris_surface_change mris_svm_classify mris_svm_train mris_twoclass mris_w_to_curv mri_threshold mri_apply_autoencoder mri_train_autoencoder mri_transform_to_COR mri_update_gca mri_vol2roi mri_wmfilter mri_xcorr oct_rf_train oct_train rbftest regdat2xfm repair_siemens_file spherical_st")

set (fake_install_script "${CMAKE_SOURCE_DIR}/attic/ATTIC_POSTINSTALL/add_dummy_binary.sh")

# If ./attic not built, then create a script for each binary during the install
# (and not build) pass, which echoes that the command is no longer built by default.
# Did not want to modify individual CMakeFiles under ./attic to  generate this script
# during the build pass when no code is compiled there. Need to pass entire command
# line as string, hence "bash -c", else execute_process will expect a subdir and fail.

if(NOT BUILD_ATTIC)
   install(CODE "
     message(STATUS \"Fake install of ./attic binaries as dummy script placeholders\")
     execute_process(COMMAND bash -c \"${fake_install_script} ${attic_binaries}\" RESULT_VARIABLE retcode)
     if(NOT \${retcode} STREQUAL 0)
       # message(FATAL_ERROR \"Install of attic placeholder scripts failed\")
       message(SEND_ERROR \"Install of attic placeholder scripts failed\")
     endif()"
   )
endif()

