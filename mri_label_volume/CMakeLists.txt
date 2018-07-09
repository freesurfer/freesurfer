project(mri_label_volume)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_label_volume mri_label_volume.c)
target_link_libraries(mri_label_volume utils)

install(TARGETS mri_label_volume DESTINATION bin)
