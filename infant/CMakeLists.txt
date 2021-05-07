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

install_pyscript(infant_recon_py)

install(FILES freesurfer_pipeline.py DESTINATION python/packages)
install(FILES CNYBCH.yaml DESTINATION average)

install(FILES license.niftyreg.txt DESTINATION docs)

install_append_help(infant_recon_all infant_recon_all.help.xml bin)

install_symlinks(all_c_robtemplate_affine_mean.05.01.2018.nii.gz TYPE files DESTINATION average)

install_pyscript(mri_label_fusion)
add_subdirectory(labelfusion)

# install external niftyreg binaries
if(INSTALL_NIFTIREG_MARTINOS)
if(APPLE)
  set(NIFTY_REG_DIR "/autofs/cluster/freesurfer/build/misc/infant/niftyreg-mac")
else()
  set(NIFTY_REG_DIR "/autofs/cluster/freesurfer/build/misc/infant/niftyreg-linux")
endif()
foreach(CMD reg_resample reg_f3d reg_aladin)
  install(PROGRAMS ${NIFTY_REG_DIR}/${CMD} DESTINATION bin)
endforeach()
endif()

# install externally-stored skullstripping models
foreach(MODEL ax_sscnn.h5 cor_sscnn.h5 sag_sscnn.h5)
  install(PROGRAMS /autofs/cluster/freesurfer/build/misc/infant/sscnn_skullstripping/${MODEL} DESTINATION average/sscnn_skullstripping)
endforeach()

# install any package requirements
if(INSTALL_PYTHON_DEPENDENCIES)
  set(PY_PKGS "numpy==1.16 tensorflow==1.5.0 keras==2.2.2 h5py==2.10.0 sklearn scikit-image tables psutil nibabel")
  if(NOT DISTRIBUTE_FSPYTHON)
    set(PKG_TARGET "--target=${CMAKE_INSTALL_PREFIX}/python/packages")
  endif()
  install(CODE "
    message(STATUS \"Installing infant-freesurfer python package dependencies\")
    execute_process(COMMAND bash -c \"${CMAKE_INSTALL_PREFIX}/python/bin/python3 -m pip install ${PKG_TARGET} ${PY_PKGS}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install infant-freesurfer python dependencies\")
    endif()
    execute_process(COMMAND bash -c \"find ${CMAKE_INSTALL_PREFIX}/python/bin ! -name 'python3*' -type f -exec rm -f {} +\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not remove pip-installed scripts\")
    endif()"
  )
endif()
