project(mris_estimate_wm)

install_pyscript(mris_estimate_wm)

install_symlinks(TYPE files DESTINATION models/topofit topofit/bary.npz)
install_symlinks(TYPE files DESTINATION models/topofit topofit/ico.npz)
install_symlinks(TYPE files DESTINATION models/topofit topofit/lh.white.ctx.average.6)
install_symlinks(TYPE files DESTINATION models/topofit topofit/mapping.npz)
install_symlinks(TYPE files DESTINATION models/topofit topofit/topofit.lh.1.pt)
install_symlinks(TYPE files DESTINATION models/topofit topofit/topofit.rh.1.pt)
