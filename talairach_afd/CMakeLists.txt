project(talairach_afd)

include_directories(${FS_INCLUDE_DIRS})

add_executable(talairach_afd talairach_afd.cpp)
add_help(talairach_afd talairach_afd.help.xml)
target_link_libraries(talairach_afd utils)

add_test_script(NAME talairach_afd_test SCRIPT test.sh)

install(TARGETS talairach_afd DESTINATION bin)
