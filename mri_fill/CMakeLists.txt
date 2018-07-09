project(mri_fill)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fill mri_fill.c)
add_help(mri_fill mri_fill.help.xml)
target_link_libraries(mri_fill utils)

install(TARGETS mri_fill DESTINATION bin)
