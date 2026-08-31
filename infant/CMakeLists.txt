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

if (FSPYTHON_INSTALL_TREE)
   set(infant_dest "${FSPYTHON_INSTALL_PREFIX}/python/packages")
   set(pydir "${FSPYTHON_INSTALL_PREFIX}/python")
else()
   set(infant_dest "python/packages")
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

if (FSPYTHON_INSTALL_TREE)
   install_pyscript_fspython_tree(infant_recon_all)
   install_pyscript_fspython_tree(mri_label_fusion)
   install(FILES freesurfer_pipeline.py DESTINATION ${infant_dest})
else()
   install_pyscript(infant_recon_all)
   install_pyscript(mri_label_fusion)
   install(FILES freesurfer_pipeline.py DESTINATION ${infant_dest})
endif()

install(FILES CNYBCH.yaml DESTINATION average)

install_symlinks(niftyreg/license.niftyreg.txt TYPE files DESTINATION docs)

#install_append_help(infant_recon_all infant_recon_all.help.xml bin)

install_symlinks(all_c_robtemplate_affine_mean.05.01.2018.nii.gz TYPE files DESTINATION average)

add_subdirectory(labelfusion)

# Entrypoint for containers
install(PROGRAMS docker/infant-container-entrypoint.bash DESTINATION bin)
install(PROGRAMS docker/infant-container-entrypoint-aws.bash DESTINATION bin)

# install external niftyreg binaries
if(APPLE)
  if(APPLE_ARM64)
    set(NIFTY_REG_DIR niftyreg/mac_arm64)
  else()
    set(NIFTY_REG_DIR niftyreg/mac_x86)
  endif()
elseif("${HOST_OS}" STREQUAL "Ubuntu20")
  set(NIFTY_REG_DIR niftyreg/ubuntu_x86)
elseif("${HOST_OS}" STREQUAL "Ubuntu22" OR "${HOST_OS}" STREQUAL "Ubuntu24")
  set(NIFTY_REG_DIR niftyreg/ubuntu22_x86)
else()
  set(NIFTY_REG_DIR niftyreg/linux_x86)
endif()

foreach(CMD reg_resample reg_f3d reg_aladin)
  install_symlinks(${NIFTY_REG_DIR}/${CMD} TYPE PROGRAMS DESTINATION bin)
endforeach()

# install the model file for skullstripping
install_symlinks(infant_synthstrip_01012025.pt TYPE files DESTINATION average/synthstrip_skullstripping)

# # add test script
# add_test_script(NAME infantfs_test SCRIPT test.sh)
