project(recon_all_clinical)

if (FSPYTHON_INSTALL_TREE)
   install_symlinks_fspython_tree(python/mri_synth_surf.py TYPE files DESTINATION python/scripts)
   install_symlinks_fspython_tree(python/norm_synthSR.py TYPE files DESTINATION python/scripts)
   install_symlinks_fspython_tree(python/refine_synthSR.py TYPE files DESTINATION python/scripts)
   install_symlinks(TYPE files DESTINATION ../fspython/models synthsurf_v10_230420.h5)
else()
   install_symlinks(python/mri_synth_surf.py TYPE files DESTINATION python/scripts)
   install_symlinks(python/norm_synthSR.py TYPE files DESTINATION python/scripts)
   install_symlinks(python/refine_synthSR.py TYPE files DESTINATION python/scripts)
   install_symlinks(TYPE files DESTINATION models synthsurf_v10_230420.h5)
endif()

install_configured(recon-all-clinical.sh DESTINATION bin)

