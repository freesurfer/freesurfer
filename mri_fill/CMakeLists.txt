project(mri_fill)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fill mri_fill.cpp)
add_help(mri_fill mri_fill.help.xml)
target_link_libraries(mri_fill utils)

install(TARGETS mri_fill DESTINATION bin)
