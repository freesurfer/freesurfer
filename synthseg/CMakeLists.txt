project(synthseg)

install_pyscript(synthseg)

install(FILES synthseg_1.0.h5 DESTINATION models)
install(FILES synthseg_segmentation_labels.npy DESTINATION models)
install(FILES synthseg_segmentation_names.npy DESTINATION models)
install(FILES synthseg_topological_classes.npy DESTINATION models)