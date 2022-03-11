project(mri_synthseg)

install_pyscript(mri_synthseg)

install(FILES synthseg_1.0.h5 DESTINATION models)
install_symlinks(TYPE files DESTINATION models synthseg_robust_1.0.h5)
install(FILES synthseg_segmentation_labels.npy DESTINATION models)
install(FILES synthseg_group_labels.npy DESTINATION models)
install(FILES synthseg_segmentation_names.npy DESTINATION models)
install(FILES synthseg_topological_classes.npy DESTINATION models)
