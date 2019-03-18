project(mri_label_vals)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_label_vals mri_label_vals.cpp)
target_link_libraries(mri_label_vals utils)

install(TARGETS mri_label_vals DESTINATION bin)
