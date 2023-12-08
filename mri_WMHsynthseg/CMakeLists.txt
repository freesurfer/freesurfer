project(mri_WMHsynthseg )

install_configured(mri_WMHsynthseg DESTINATION bin)

install(DIRECTORY WMHSynthSeg DESTINATION python/packages)
