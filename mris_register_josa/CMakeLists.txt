project(mris_register_josa)

if (FSPYTHON_INSTALL_TREE)
   set(mris_register_josa_dest "${FSPYTHON_INSTALL_PREFIX}/models")
else()
   set(mris_register_josa_dest "models")
endif()

install(FILES mris_register_josa_20241121_lh.h5 DESTINATION ${mris_register_josa_dest})
install(FILES mris_register_josa_20241121_rh.h5 DESTINATION ${mris_register_josa_dest})

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mris_register_josa)
else()
   install_pyscript(mris_register_josa)
endif()

