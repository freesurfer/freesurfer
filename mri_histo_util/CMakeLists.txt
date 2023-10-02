project(mri_histo_atlas_segment)

install_configured(mri_histo_atlas_segment DESTINATION bin)

install(DIRECTORY ERC_bayesian_segmentation DESTINATION python/packages)

