project(mri_WMHsynthseg )

install_configured(mri_WMHsynthseg DESTINATION bin)

if (FSPYTHON_INSTALL_TREE)
   install(DIRECTORY WMHSynthSeg DESTINATION ../fspython/python/packages)
else()
   install(DIRECTORY WMHSynthSeg DESTINATION python/packages)
endif()

