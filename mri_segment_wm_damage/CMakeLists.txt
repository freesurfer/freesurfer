project(mri_segment_wm_damage)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_segment_wm_damage mri_segment_wm_damage.c)
target_link_libraries(mri_segment_wm_damage utils)

install(TARGETS mri_segment_wm_damage DESTINATION bin)
