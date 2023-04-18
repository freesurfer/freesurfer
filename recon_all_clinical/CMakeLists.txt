project(recon_all_clinical)

install_symlinks(TYPE files DESTINATION python mri_synt_surf.py)
install_symlinks(TYPE files DESTINATION python mri_synthSR.py)
install_symlinks(TYPE files DESTINATION python refine_synthSR.py)


install_pyscript(mri_synt_surf)
install_pyscript(mri_synthSR)
install_pyscript(refine_synthSR)

install(FILES mri_model.h5)