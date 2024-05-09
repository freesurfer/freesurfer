project(exvivo)

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_exvivo_strip)
   install_pyscript_fspython_tree(mri_exvivo_norm)
   install_symlinks(TYPE files DESTINATION ../fspython/models exvivo.strip.lh.h5 exvivo.norm.lh.h5)
else()
   install_pyscript(mri_exvivo_strip)
   install_pyscript(mri_exvivo_norm)
   install_symlinks(TYPE files DESTINATION models exvivo.strip.lh.h5 exvivo.norm.lh.h5)
endif()

