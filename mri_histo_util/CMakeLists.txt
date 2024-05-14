project(mri_histo_atlas_segment)

install_configured(mri_histo_atlas_segment DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install(DIRECTORY ERC_bayesian_segmentation DESTINATION ../fspython/python/packages)
else()
   install(DIRECTORY ERC_bayesian_segmentation DESTINATION python/packages)
endif()

