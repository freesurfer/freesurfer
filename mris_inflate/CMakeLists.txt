project(mris_inflate)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_inflate mris_inflate.c)
add_help(mris_inflate mris_inflate.help.xml)
target_link_libraries(mris_inflate utils)

install(TARGETS mris_inflate DESTINATION bin)
