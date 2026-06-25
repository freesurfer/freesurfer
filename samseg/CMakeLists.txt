project(samseg)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${FSPYTHON_INSTALL_PREFIX}/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()


if(NOT INTEGRATE_SAMSEG)
   message(STATUS "No install of run_samseg run_samseg_long sbtiv gems_compute_atlas_probs for INTEGRATE_SAMSEG=OFF")
   # if (FSPYTHON_INSTALL_TREE)
   #    install_pyscript_fspython_tree(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   # else()
   #   install_pyscript(run_samseg run_samseg_long sbtiv gems_compute_atlas_probs)
   # endif()
# else()
   ## the case for -DINTEGRATE_SAMSEG=ON is handled by ./samseg_scripts/CMakeLists.txt
   # install_samseg_cli()   
endif()


install_configured(samseg samseg-long DESTINATION bin)

## should work on intel and silicon Mac (12/29/2025)
# if(NOT APPLE_ARM64)
## turn off for darwin_x86_64 (05/14/2026)
if(NOT APPLE)
   add_test_script(NAME samseg_test SCRIPT test.sh)
elseif(APPLE_ARM64)
   add_test_script(NAME samseg_test SCRIPT test.sh)
endif()
# endif()

