project(mri_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segment mri_segment.cpp)
add_help(mri_segment mri_segment.help.xml)
target_link_libraries(mri_segment utils)

add_test_script(NAME mri_segment_test SCRIPT test.sh)

install(TARGETS mri_segment DESTINATION bin)
