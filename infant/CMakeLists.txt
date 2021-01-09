project(infant)

install(PROGRAMS
  atlassubjects.csh
  check_infant_recon_success.csh
  create_pial_surfaces_mprage_subject.csh
  create_wm_surfaces_mprage_subject.csh
  find_cerebral_wm_labels.csh
  find_neighbors_byAge.csh
  find_neighbors_byAge_extended.csh
  find_neighbors_byMI.csh
  find_neighbors_byMI_extended.csh
  run_label_fusion.csh
  set_babydev_packages.csh
  set_infant_recon_params.csh
  DESTINATION bin
)

install_append_help(infant_recon_all infant_recon_all.help.xml bin)

install_symlinks(all_c_robtemplate_affine_mean.05.01.2018.nii.gz TYPE files DESTINATION average)

install_pyscript(mri_label_fusion)
add_subdirectory(labelfusion)
