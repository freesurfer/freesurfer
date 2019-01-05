project(mris_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_register mris_register.c)
add_help(mris_register mris_register.help.xml)
target_link_libraries(mris_register utils)

add_test_script(NAME mris_register_test SCRIPT test_mris_register DEPENDS mris_register)
configure_file(testdata.tar.gz . COPYONLY)

install(TARGETS mris_register DESTINATION bin)
