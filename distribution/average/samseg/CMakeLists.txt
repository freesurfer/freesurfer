# samseg atlases

install_directories(

    DESTINATION average/samseg

    20Subjects_smoothing2_down2_smoothingForAffine2

    20Subjects_smoothing2_down2_smoothingForAffine2_lesion

    20Subjects_smoothing2_down2_smoothingForAffine2_lesion_wm_prior

)

 

install_symlinks(
  TYPE files
  DESTINATION average/samseg
  samseg.talairach.lta.mrf.gca
  samseg.talairach.m3z.mrf.gca
)
