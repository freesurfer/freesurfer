project(mri_label_volume)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label_volume mri_label_volume.cpp)
target_link_libraries(mri_label_volume utils)

install(TARGETS mri_label_volume DESTINATION bin)
