project(mri_concatenate_lta)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_concatenate_lta mri_concatenate_lta.cpp)
target_link_libraries(mri_concatenate_lta utils)

add_test_script(NAME mri_concatenate_lta_test SCRIPT test.sh)

install(TARGETS mri_concatenate_lta DESTINATION bin)