project(recon_all_clinical)

install_symlinks(python/mri_synth_surf.py TYPE files DESTINATION python/scripts)
install_symlinks(python/norm_synthSR.py TYPE files DESTINATION python/scripts)
install_symlinks(python/refine_synthSR.py TYPE files DESTINATION python/scripts)

install_configured(recon-all-clinical.sh DESTINATION bin)

install_symlinks(TYPE files DESTINATION models synthsurf_v10_230420.h5)
