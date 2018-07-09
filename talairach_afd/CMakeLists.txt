project(talairach_afd)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(talairach_afd talairach_afd.c)
add_help(talairach_afd talairach_afd.help.xml)
target_link_libraries(talairach_afd utils)

install(TARGETS talairach_afd DESTINATION bin)
