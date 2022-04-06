project(mri_concatenate_gcam)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_concatenate_gcam mri_concatenate_gcam.cpp)
add_help(mri_concatenate_gcam mri_concatenate_gcam.help.xml)
target_link_libraries(mri_concatenate_gcam utils)

add_test_script(NAME mri_concatenate_gcam_test SCRIPT test.sh)

install(TARGETS mri_concatenate_gcam DESTINATION bin)
