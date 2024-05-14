project(mri_synthseg)

if (FSPYTHON_INSTALL_TREE)
   set(mri_synthsr_dest "\.\./fspython/models")
else()
   set(mri_synthsr_dest "models")
endif()

install(FILES synthseg_1.0.h5 DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_segmentation_labels.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_labels_2.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_names.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_names_2.0.npy DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_denoiser_labels_2.0.npy DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_parcellation_labels.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_parcellation_names.npy DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_qc_labels.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_qc_labels_2.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_qc_names.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_qc_names_2.0.npy DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_topological_classes.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_topological_classes_2.0.npy DESTINATION ${mri_synthsr_dest})

install(FILES synthseg_segmentation_labels_photo_both_1.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_labels_photo_single_1.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_names_photo_both_1.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_segmentation_names_photo_single_1.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_topological_classes_photo_both_1.0.npy DESTINATION ${mri_synthsr_dest})
install(FILES synthseg_topological_classes_photo_single_1.0.npy DESTINATION ${mri_synthsr_dest})

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_synthseg)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_robust_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_parc_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_qc_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_photo_both_1.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_photo_single_1.0.h5)
else()
   install_pyscript(mri_synthseg)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_2.0.h5)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_robust_2.0.h5)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_parc_2.0.h5)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_qc_2.0.h5)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_photo_both_1.0.h5)
   install_symlinks(TYPE files DESTINATION ${mri_synthsr_dest} synthseg_photo_single_1.0.h5)
endif()

