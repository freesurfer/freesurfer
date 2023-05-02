project(recon_all_clinical)

install_symlinks(TYPE files DESTINATION python mri_synth_surf.py)
install_symlinks(TYPE files DESTINATION python norm_synthSR.py)
install_symlinks(TYPE files DESTINATION python refine_synthSR.py)

install_configured(recon-all-clinical.sh DESTINATION bin)

install_symlinks(TYPE files DESTINATION models synthsurf_v10_230420.h5)
