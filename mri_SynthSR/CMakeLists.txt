project(mri_SynthSR)

install_pyscript(mri_SynthSR)
install_pyscript(mri_SynthSR_hyperfine)

install(FILES SynthSR_v10_210712.h5 DESTINATION models)
install(FILES SynthSR_v10_210712_hyperfine.h5 DESTINATION models)

