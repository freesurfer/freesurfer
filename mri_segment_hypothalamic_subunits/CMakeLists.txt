project(mri_segment_hypothalamic_subunits)

if (FSPYTHON_INSTALL_TREE)
   set(mri_synthsr_dest "\.\./fspython/models")
else()
   set(mri_synthsr_dest "models")
endif()

install(FILES hypothalamic_subunits.h5 DESTINATION ${mri_synthsr_dest})

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_segment_hypothalamic_subunits)
else()
   install_pyscript(mri_segment_hypothalamic_subunits)
endif()

