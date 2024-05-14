project(samseg)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${CMAKE_INSTALL_PREFIX}/\.\./fspython/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

function(install_samsegscript)
  # ??? todo: retrieve samseg directory using get_package_dir(samseg samsegdir) ???
  # set SAMSEGCLI_DIR, this is hard-coded for now
  if(NOT "${USER_SAMSEG_PATH}" STREQUAL "")
    set(SAMSEGCLI_DIR ${USER_SAMSEG_PATH})
  else()
    set(SAMSEGCLI_DIR ${pydir}/lib/python3.8/site-packages)
    if(NOT DISTRIBUTE_FSPYTHON)
      set(SAMSEGCLI_DIR ${pydir}/packages)
    endif()
  endif()

  # config the scripts
  foreach(SCRIPT ${ARGN})
    # config fspython wrapper script
    message(STATUS "install processing samseg script ${SCRIPT}")
    if (FSPYTHON_INSTALL_TREE)
       config_fspythonwrapper_fspython_tree(${SCRIPT})
    else()
       config_fspythonwrapper(${SCRIPT})
    endif()

    # remove pip-installed cli wrappers
    if(DISTRIBUTE_FSPYTHON)
      remove_files(${pydir}/bin/${SCRIPT})
    else()
      remove_files(${pydir}/packages/bin/${SCRIPT})
    endif()

    # create links in ${pydir}/scripts/ to real python scripts in installed package
    symlink(${SAMSEGCLI_DIR}/samseg/cli/${SCRIPT}.py ${pydir}/scripts/${SCRIPT})
  endforeach()
endfunction()

install_configured(samseg samseg-long DESTINATION bin)
if(NOT INTEGRATE_SAMSEG)
   if (FSPYTHON_INSTALL_TREE)
      install_pyscript_fspython_tree(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   else()
      install_pyscript(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
  endif()
  add_test_script(NAME samseg SCRIPT test.sh)
else()
  install_samsegscript(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs segment_subregions prepareAtlasDirectory merge_add_mesh_alphas computeTissueConcentrations)
  add_test(NAME samseg COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test.sh)
endif()

