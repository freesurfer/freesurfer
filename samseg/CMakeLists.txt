project(samseg)

function(install_samsegscript)
  # ??? todo: retrieve samseg directory using get_package_dir(samseg samsegdir) ???
  # set SAMSEGCLI_DIR, this is hard-coded for now
  if(NOT "${USER_SAMSEG_PATH}" STREQUAL "")
    set(SAMSEGCLI_DIR ${USER_SAMSEG_PATH})
  else()
    set(SAMSEGCLI_DIR ${CMAKE_INSTALL_PREFIX}/python/lib/python3.8/site-packages)
    if(NOT DISTRIBUTE_FSPYTHON)
      set(SAMSEGCLI_DIR ${CMAKE_INSTALL_PREFIX}/python/packages)
    endif()
  endif()

  # config the scripts
  foreach(SCRIPT ${ARGN})
    # config fspython wrapper script
    config_fspythonwrapper(${SCRIPT})

    # remove pip-installed cli wrappers
    if(DISTRIBUTE_FSPYTHON)
      remove_files(${CMAKE_INSTALL_PREFIX}/python/bin/${SCRIPT})
    else()
      remove_files(${CMAKE_INSTALL_PREFIX}/python/packages/bin/${SCRIPT})
    endif()

    # create links in ${CMAKE_INSTALL_PREFIX}/python/scripts/ to real python scripts in installed package
    symlink(${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py ${CMAKE_INSTALL_PREFIX}/python/scripts/${SCRIPT})
  endforeach()
endfunction()


install_configured(samseg samseg-long DESTINATION bin)
if(NOT INTEGRATE_SAMSEG)
  install_pyscript(
    run_samseg
    run_samseg_long
    sbtiv
    gems_compute_atlas_probs
  )

  add_test_script(NAME samseg SCRIPT test.sh)
else()
  install_samsegscript(
    run_samseg
    run_samseg_long
    sbtiv
    gems_compute_atlas_probs
    segment_subregions
    prepareAtlasDirectory
    merge_add_mesh_alphas
    computeTissueConcentrations    
  )
  
  add_test(NAME samseg COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test.sh)
endif()





