project(mri_strip_subject_info)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_strip_subject_info mri_strip_subject_info.cpp)
target_link_libraries(mri_strip_subject_info utils)

install(TARGETS mri_strip_subject_info DESTINATION bin)
