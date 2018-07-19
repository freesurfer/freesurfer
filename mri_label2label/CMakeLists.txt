project(mri_label2label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_label2label mri_label2label.c)
target_link_libraries(mri_label2label utils)
install_wrapped(TARGETS mri_label2label DESTINATION bin)
