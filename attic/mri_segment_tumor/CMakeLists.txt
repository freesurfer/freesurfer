project(mri_segment_tumor)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segment_tumor mri_segment_tumor.cpp)
target_link_libraries(mri_segment_tumor utils)

install(TARGETS mri_segment_tumor DESTINATION bin)
