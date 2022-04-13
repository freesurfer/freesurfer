project(mri_label2vol)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label2vol mri_label2vol.cpp)
target_link_libraries(mri_label2vol utils)

add_test_script(NAME mri_label2vol_test SCRIPT test.sh)

install(TARGETS mri_label2vol DESTINATION bin)
