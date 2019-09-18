project(mri_segment_wm_damage)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segment_wm_damage mri_segment_wm_damage.cpp)
target_link_libraries(mri_segment_wm_damage utils)

install(TARGETS mri_segment_wm_damage DESTINATION bin)
