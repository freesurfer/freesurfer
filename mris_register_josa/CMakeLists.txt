project(mris_register_josa)

install_pyscript(mris_register_josa)

install(FILES mris_register_josa_lh_20231025.h5 DESTINATION models)

install_directories(spheremorph DESTINATION python/packages)