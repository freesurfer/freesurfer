project(mri_synthsr)

install_pyscript(mri_synthsr)
install_pyscript(mri_synthsr_hyperfine)

install(FILES synthsr_v10_210712.h5 DESTINATION models)
install(FILES synthsr_v10_210712_hyperfine.h5 DESTINATION models)