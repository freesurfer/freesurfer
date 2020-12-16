# samseg atlases

install_tarball(20Subjects_smoothing2_down2_smoothingForAffine2.tar.gz DESTINATION average/samseg)
install_tarball(20Subjects_smoothing2_down2_smoothingForAffine2_lesion.tar.gz DESTINATION average/samseg)
install_tarball(20Subjects_smoothing2_down2_smoothingForAffine2_lesion_wm_prior.tar.gz DESTINATION average/samseg)

install_symlinks(
  TYPE files
  DESTINATION average/samseg
  samseg.talairach.lta.mrf.gca
  samseg.talairach.m3z.mrf.gca
)
