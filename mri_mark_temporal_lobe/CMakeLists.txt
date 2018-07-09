project(mri_mark_temporal_lobe)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_mark_temporal_lobe mri_mark_temporal_lobe.c)
target_link_libraries(mri_mark_temporal_lobe utils)

install(TARGETS mri_mark_temporal_lobe DESTINATION bin)
