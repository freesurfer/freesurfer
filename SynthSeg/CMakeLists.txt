project(mri_segment_hypothalamic_subunits)

install_pyscript(SynthSeg)

install(FILES SynthSeg_1.0.h5 DESTINATION models)
install(FILES SynthSeg_segmentation_labels.npy DESTINATION models)
install(FILES SynthSeg_segmentation_names.npy DESTINATION models)
install(FILES SynthSeg_topological_classes.npy DESTINATION models)