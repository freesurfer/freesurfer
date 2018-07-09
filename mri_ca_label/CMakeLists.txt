project(mri_ca_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ca_label mri_ca_label.c)
add_help(mri_ca_label mri_ca_label.help.xml)
target_link_libraries(mri_ca_label utils)

install(TARGETS mri_ca_label DESTINATION bin)
