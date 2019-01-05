project(mri_ca_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ca_label mri_ca_label.c)
add_help(mri_ca_label mri_ca_label.help.xml)
target_link_libraries(mri_ca_label utils)

add_test_script(NAME mri_ca_label_test SCRIPT test_mri_ca_label DEPENDS mri_ca_label)
configure_file(testdata.tar.gz . COPYONLY)

install(TARGETS mri_ca_label DESTINATION bin)
