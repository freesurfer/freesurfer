project(mri_convert_mdh)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_convert_mdh mri_convert_mdh.c)
target_link_libraries(mri_convert_mdh utils)

install(TARGETS mri_convert_mdh DESTINATION bin)
