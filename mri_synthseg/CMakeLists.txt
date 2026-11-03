project(mri_synthseg)

if (FSPYTHON_INSTALL_TREE)
   set(mri_synthsr_dest "${FSPYTHON_INSTALL_PREFIX}/models")
else()
   set(mri_synthsr_dest "models")
endif()

# use a soft link for mri_synthseg to be included in the fspylegacy as synthseg.py
ADD_CUSTOM_TARGET(mri_synthseg ALL COMMAND ${CMAKE_COMMAND} -E create_symlink \"${CMAKE_SOURCE_DIR}/mri_synthseg/mri_synthseg\" \"${CMAKE_SOURCE_DIR}/python/fspylegacy/synthseg.py\")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR}/mri_synthseg APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${CMAKE_SOURCE_DIR}/python/fspylegacy/synthseg.py)

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
   install_modelpredict_pyscript_fspython_tree(synthseg)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_robust_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_parc_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_qc_2.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_photo_both_1.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthseg_photo_single_1.0.h5)
else()
   install_modelpredict_pyscript(synthseg)
   install_symlinks(TYPE files DESTINATION models synthseg_2.0.h5)
   install_symlinks(TYPE files DESTINATION models synthseg_robust_2.0.h5)
   install_symlinks(TYPE files DESTINATION models synthseg_parc_2.0.h5)
   install_symlinks(TYPE files DESTINATION models synthseg_qc_2.0.h5)
   install_symlinks(TYPE files DESTINATION models synthseg_photo_both_1.0.h5)
   install_symlinks(TYPE files DESTINATION models synthseg_photo_single_1.0.h5)
endif()

