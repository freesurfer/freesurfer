project(mri_make_template)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_make_template mri_make_template.c)
target_link_libraries(mri_make_template utils)

install(TARGETS mri_make_template DESTINATION bin)
