project(mri_extract_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_extract_label mri_extract_label.cpp)
target_link_libraries(mri_extract_label utils)

install(TARGETS mri_extract_label DESTINATION bin)
