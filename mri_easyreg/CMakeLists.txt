project(mri_easyreg)

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_easyreg mri_easywarp)
   install_symlinks_fspython_tree(TYPE files DESTINATION models easyreg_v10_230103.h5)
else() 
   install_pyscript(mri_easyreg mri_easywarp)
   install_symlinks(TYPE files DESTINATION models easyreg_v10_230103.h5)
endif()

