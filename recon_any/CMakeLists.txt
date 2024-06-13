project(recon_any)

install_configured(run_recon-any DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
    install(DIRECTORY recon-any DESTINATION ../fspython/python/packages)
else()
    install(DIRECTORY recon-any DESTINATION python/packages)
endif()