project(mris_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_register mris_register.cpp)
add_help(mris_register mris_register.help.xml)
target_link_libraries(mris_register utils)

add_test_script(NAME mris_register_test SCRIPT test.sh)

install(TARGETS mris_register DESTINATION bin)
