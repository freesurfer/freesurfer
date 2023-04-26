project(mri_synthseg)

install_pyscript(mri_synthseg)

install(FILES synthseg_1.0.h5 DESTINATION models)
install_symlinks(TYPE files DESTINATION models synthseg_2.0.h5)
install_symlinks(TYPE files DESTINATION models synthseg_robust_2.0.h5)
install_symlinks(TYPE files DESTINATION models synthseg_parc_2.0.h5)
install_symlinks(TYPE files DESTINATION models synthseg_qc_2.0.h5)

install(FILES synthseg_segmentation_labels.npy DESTINATION models)
install(FILES synthseg_segmentation_labels_2.0.npy DESTINATION models)
install(FILES synthseg_segmentation_names.npy DESTINATION models)
install(FILES synthseg_segmentation_names_2.0.npy DESTINATION models)

install(FILES synthseg_denoiser_labels_2.0.npy DESTINATION models)

install(FILES synthseg_parcellation_labels.npy DESTINATION models)
install(FILES synthseg_parcellation_names.npy DESTINATION models)

install(FILES synthseg_qc_labels.npy DESTINATION models)
install(FILES synthseg_qc_labels_2.0.npy DESTINATION models)
install(FILES synthseg_qc_names.npy DESTINATION models)
install(FILES synthseg_qc_names_2.0.npy DESTINATION models)

install(FILES synthseg_topological_classes.npy DESTINATION models)
install(FILES synthseg_topological_classes_2.0.npy DESTINATION models)


install_symlinks(TYPE files DESTINATION models synthseg_photo_both_1.0.h5)
install_symlinks(TYPE files DESTINATION models synthseg_photo_single_1.0.h5)
install(FILES synthseg_segmentation_labels_photo_both_1.0.npy DESTINATION models)
install(FILES synthseg_segmentation_labels_photo_single_1.0.npy DESTINATION models)
install(FILES synthseg_segmentation_names_photo_both_1.0.npy DESTINATION models)
install(FILES synthseg_segmentation_names_photo_single_1.0.npy DESTINATION models)
install(FILES synthseg_topological_classes_photo_both_1.0.npy DESTINATION models)
install(FILES synthseg_topological_classes_photo_single_1.0.npy DESTINATION models)