project(mri_label_vals)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_label_vals mri_label_vals.c)
target_link_libraries(mri_label_vals utils)

install(TARGETS mri_label_vals DESTINATION bin)
