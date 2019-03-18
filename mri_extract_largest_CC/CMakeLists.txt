project(mri_extract_largest_CC)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_extract_largest_CC mri_extract_largest_CC.cpp)
target_link_libraries(mri_extract_largest_CC utils)

install(TARGETS mri_extract_largest_CC DESTINATION bin)
