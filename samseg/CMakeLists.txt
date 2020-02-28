project(samseg)

install(PROGRAMS samseg samseg-long DESTINATION bin)
install_pyscript(run_samseg run_samseg_long sbtiv)
