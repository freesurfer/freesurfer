project(mri_ca_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ca_register mri_ca_register.c)
add_help(mri_ca_register mri_ca_register.help.xml)
target_link_libraries(mri_ca_register utils)

install(TARGETS mri_ca_register DESTINATION bin)
