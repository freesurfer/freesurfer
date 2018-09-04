project(lta_convert)

include_directories(${FS_INCLUDE_DIRS})

add_executable(lta_convert lta_convert.cpp)
add_help(lta_convert lta_convert.help.xml)
target_link_libraries(lta_convert utils)

install(TARGETS lta_convert DESTINATION bin)
