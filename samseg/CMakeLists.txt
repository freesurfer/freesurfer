project(samseg)

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
  integrate_samseg(
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





