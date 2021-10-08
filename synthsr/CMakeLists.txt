project(synthsr)

install_pyscript(synthsr)
install_pyscript(synthsr_hyperfine)

install(FILES synthsr_v10_210712.h5 DESTINATION models)
install(FILES synthsr_v10_210712_hyperfine.h5 DESTINATION models)