project(mris_volume)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_volume mris_volume.c)
target_link_libraries(mris_volume utils)
install(TARGETS mris_volume DESTINATION bin)

add_executable(mris_wm_volume mris_wm_volume.c)
target_link_libraries(mris_wm_volume utils)
install(TARGETS mris_wm_volume DESTINATION bin)

