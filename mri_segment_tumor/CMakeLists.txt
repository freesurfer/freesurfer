project(mri_segment_tumor)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_segment_tumor mri_segment_tumor.c)
target_link_libraries(mri_segment_tumor utils)

install(TARGETS mri_segment_tumor DESTINATION bin)
