project(mri_mark_temporal_lobe)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mark_temporal_lobe mri_mark_temporal_lobe.cpp)
target_link_libraries(mri_mark_temporal_lobe utils)

install(TARGETS mri_mark_temporal_lobe DESTINATION bin)
