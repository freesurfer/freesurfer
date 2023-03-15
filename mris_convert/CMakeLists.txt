project(mris_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_convert mris_convert.cpp)
add_help(mris_convert mris_convert.help.xml)
target_link_libraries(mris_convert utils)

add_test_script(NAME mris_convert_test SCRIPT test.sh)

install(TARGETS mris_convert DESTINATION bin)
