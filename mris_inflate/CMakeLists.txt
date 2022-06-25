project(mris_inflate)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_inflate mris_inflate.cpp)
add_help(mris_inflate mris_inflate.help.xml)
target_link_libraries(mris_inflate utils)

add_test_script(NAME mris_inflate_posix_test SCRIPT test.sh)

install(TARGETS mris_inflate DESTINATION bin)
