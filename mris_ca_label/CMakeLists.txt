project(mris_ca_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_ca_label mris_ca_label.cpp)
add_help(mris_ca_label mris_ca_label.help.xml)
target_link_libraries(mris_ca_label utils)

add_test_script(NAME mris_ca_label_test SCRIPT test.sh)

install(TARGETS mris_ca_label DESTINATION bin)
