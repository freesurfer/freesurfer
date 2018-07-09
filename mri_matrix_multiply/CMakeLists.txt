project(mri_matrix_multiply)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_matrix_multiply mri_matrix_multiply.c)
target_link_libraries(mri_matrix_multiply utils)

install(TARGETS mri_matrix_multiply DESTINATION bin)
