project(mri_add_xform_to_header)

include_directories(${CMAKE_SOURCE_DIR}/include ${MINC_INCLUDE_DIR})

add_executable(mri_add_xform_to_header mri_add_xform_to_header.c)
add_help(mri_add_xform_to_header mri_add_xform_to_header.help.xml)
target_link_libraries(mri_add_xform_to_header utils)

install(TARGETS mri_add_xform_to_header DESTINATION bin)
