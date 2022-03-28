project(exvivo)

install_pyscript(mri_exvivo_strip)
install_pyscript(mri_exvivo_norm)

install_symlinks(
    TYPE files DESTINATION models
    exvivo.strip.lh.h5
    exvivo.norm.lh.h5
)
