project(mri_synthsr)

if (FSPYTHON_INSTALL_TREE)
   set(mri_synthsr_dest "\.\./fspython/models")
else()
   set(mri_synthsr_dest "models")
endif()

# add_custom_target(check_h5_files ALL COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/h5.sh ${CMAKE_SOURCE_DIR})
install(FILES synthsr_v10_210712.h5 DESTINATION ${mri_synthsr_dest})
install(FILES synthsr_v10_210712_hyperfine.h5 DESTINATION ${mri_synthsr_dest})

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_synthsr)
   install_pyscript_fspython_tree(mri_synthsr_hyperfine)
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthsr_v20_230130.h5 )
   install_symlinks_fspython_tree(TYPE files DESTINATION models synthsr_lowfield_v20_230130.h5)
else()
   install_pyscript(mri_synthsr)
   install_pyscript(mri_synthsr_hyperfine)
   install_symlinks(TYPE files DESTINATION models synthsr_v20_230130.h5 )
   install_symlinks(TYPE files DESTINATION models synthsr_lowfield_v20_230130.h5)
endif()

add_test_script(NAME mri_synthsr_test SCRIPT test.sh)

