project(mri_make_template)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_make_template mri_make_template.cpp)
target_link_libraries(mri_make_template utils)

install(TARGETS mri_make_template DESTINATION bin)
