project(mri_extract_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_extract_label mri_extract_label.c)
target_link_libraries(mri_extract_label utils)

install(TARGETS mri_extract_label DESTINATION bin)
