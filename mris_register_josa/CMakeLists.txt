project(mris_register_josa)

install_pyscript(mris_register_josa)

install(FILES mris_register_josa_rh_20231025.h5 DESTINATION models)
install(FILES mris_register_josa_lh_20231025.h5 DESTINATION models)

install_directories(spheremorph/spheremorph DESTINATION python/packages)