project(mri_matrix_multiply)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_matrix_multiply mri_matrix_multiply.cpp)
target_link_libraries(mri_matrix_multiply utils)

add_test_script(NAME mri_matrix_multiply_test SCRIPT test.sh)

install(TARGETS mri_matrix_multiply DESTINATION bin)
