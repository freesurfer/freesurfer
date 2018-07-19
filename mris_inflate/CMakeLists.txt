project(mris_inflate)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_inflate mris_inflate.c)
add_help(mris_inflate mris_inflate.help.xml)
target_link_libraries(mris_inflate utils)

install(TARGETS mris_inflate DESTINATION bin)
