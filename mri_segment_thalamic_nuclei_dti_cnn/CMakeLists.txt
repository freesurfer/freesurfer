project(mri_segment_thalamic_nuclei_dti_cnn)

if (FSPYTHON_INSTALL_TREE)
   set(thalseg_dest "\.\./fspython/models")
   install_pyscript_fspython_tree(mri_segment_thalamic_nuclei_dti_cnn)
   install_symlinks_fspython_tree(TYPE files DESTINATION models thalseg_1.0.h5)
   install_symlinks_fspython_tree(TYPE files DESTINATION models thalseg_1.1.h5)
else()
   set(thalseg_dest "models")
   install_pyscript(mri_segment_thalamic_nuclei_dti_cnn)
   install_symlinks(TYPE files DESTINATION models thalseg_1.0.h5)
   install_symlinks(TYPE files DESTINATION models thalseg_1.1.h5)
endif()

install(FILES thalseg_segmentation_labels_1.0.npy DESTINATION ${thalseg_dest})
install(FILES thalseg_segmentation_names_1.0.npy DESTINATION ${thalseg_dest})
install(FILES thalseg_topological_classes_1.0.npy DESTINATION ${thalseg_dest})

