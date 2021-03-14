project(mri_cvs_register)

install_append_help(mri_cvs_register mri_cvs_register.help.xml bin)
install_configured(mri_cvs_check mri_cvs_data_copy DESTINATION bin)
install(FILES mri_cvs_requiredfiles.txt id.xfm DESTINATION bin)
