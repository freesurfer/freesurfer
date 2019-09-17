project(mri_make_labels)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_make_labels mri_make_labels.cpp)
target_link_libraries(mri_make_labels utils)

install(TARGETS mri_make_labels DESTINATION bin)
