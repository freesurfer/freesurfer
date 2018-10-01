project(template)

include_directories(${FS_INCLUDE_DIRS})

add_executable(template template.c)
target_link_libraries(template utils)

install(TARGETS template DESTINATION bin)
