project(mri_segmentThalamicNucleiDTI_CNN)

install_pyscript(mri_segmentThalamicNucleiDTI_CNN)

install_symlinks(TYPE files DESTINATION models thalseg_1.0.h5)

install(FILES thalseg_segmentation_labels_1.0.npy DESTINATION models)
install(FILES thalseg_segmentation_names_1.0.npy DESTINATION models)
install(FILES thalseg_topological_classes_1.0.npy DESTINATION models)

