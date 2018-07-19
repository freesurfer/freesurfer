project(talairach_afd)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(talairach_afd talairach_afd.c)
add_help(talairach_afd talairach_afd.help.xml)
target_link_libraries(talairach_afd utils)

add_test_script(NAME talairach_afd_test SCRIPT test_tal_afd DEPENDS talairach_afd)

install(TARGETS talairach_afd DESTINATION bin)
