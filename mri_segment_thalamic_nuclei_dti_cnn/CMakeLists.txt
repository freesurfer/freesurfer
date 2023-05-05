project(mri_segment_thalamic_nuclei_dti_cnn)

install_pyscript(mri_segment_thalamic_nuclei_dti_cnn)

install_symlinks(TYPE files DESTINATION models thalseg_1.0.h5)

install(FILES thalseg_segmentation_labels_1.0.npy DESTINATION models)
install(FILES thalseg_segmentation_names_1.0.npy DESTINATION models)
install(FILES thalseg_topological_classes_1.0.npy DESTINATION models)

