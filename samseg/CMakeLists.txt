project(samseg)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)

#
add_executable(samsegmesh2surf samsegmesh2surf.cpp)
target_link_libraries(samsegmesh2surf kvlGEMSCommon utils)

install_configured(samseg samseg-long DESTINATION bin)
install_pyscript(
  run_samseg
  run_samseg_long
  sbtiv
  gems_compute_atlas_probs
)

add_test_script(NAME samseg SCRIPT test.sh)
