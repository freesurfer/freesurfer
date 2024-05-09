project(mri_3d_photo_recon)

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(mri_3d_photo_recon)
else()
   install_pyscript(mri_3d_photo_recon)
endif()

