project(mris_register_josa)

if (FSPYTHON_INSTALL_TREE)
   set(mri_synthsr_dest "\.\./fspython/models")
else()
   set(mri_synthsr_dest "models")
endif()

install(FILES mris_register_josa_20240401_rh.h5 DESTINATION ${mri_synthsr_dest})
install(FILES mris_register_josa_20240401_lh.h5 DESTINATION ${mri_synthsr_dest})

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mris_register_josa)
   install_directories_fspython_tree(spheremorph DESTINATION python/packages)
else()
   install_pyscript(mris_register_josa)
   install_directories(spheremorph DESTINATION python/packages)
endif()

