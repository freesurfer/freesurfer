project(mri_synthsr)

install_pyscript(mri_synthsr)
install_pyscript(mri_synthsr_hyperfine)

# add_custom_target(check_h5_files ALL COMMAND bash ${CMAKE_SOURCE_DIR}/cmake/h5.sh ${CMAKE_SOURCE_DIR})
install(FILES synthsr_v10_210712.h5 DESTINATION models)
install(FILES synthsr_v10_210712_hyperfine.h5 DESTINATION models)

install_symlinks(TYPE files DESTINATION models synthsr_v20_230130.h5 )
install_symlinks(TYPE files DESTINATION models synthsr_lowfield_v20_230130.h5)

add_test_script(NAME mri_synthsr_test SCRIPT test.sh)

