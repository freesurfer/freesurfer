project(mri_extract_largest_CC)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_extract_largest_CC mri_extract_largest_CC.c)
target_link_libraries(mri_extract_largest_CC utils)

install(TARGETS mri_extract_largest_CC DESTINATION bin)
