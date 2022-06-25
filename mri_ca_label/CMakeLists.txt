project(mri_ca_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_label mri_ca_label.cpp)
add_help(mri_ca_label mri_ca_label.help.xml)
target_link_libraries(mri_ca_label utils)

add_test_script(NAME mri_ca_label_test SCRIPT test.sh)

install(TARGETS mri_ca_label DESTINATION bin)
