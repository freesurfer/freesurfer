project(mri_mask)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_mask mri_mask.c)
add_help(mri_mask mri_mask.help.xml)
target_link_libraries(mri_mask utils)
install_wrapped(TARGETS mri_mask DESTINATION bin)
