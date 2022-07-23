project(mri_ca_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_register mri_ca_register.cpp)
add_help(mri_ca_register mri_ca_register.help.xml)
target_link_libraries(mri_ca_register utils)

add_test_script(NAME mri_ca_register_test SCRIPT test.sh)

install(TARGETS mri_ca_register DESTINATION bin)
