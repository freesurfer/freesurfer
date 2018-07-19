project(mri_make_labels)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_make_labels mri_make_labels.c)
target_link_libraries(mri_make_labels utils)

install(TARGETS mri_make_labels DESTINATION bin)
