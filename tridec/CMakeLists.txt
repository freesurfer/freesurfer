project(tridec)

include_directories(${FS_INCLUDE_DIRS})

add_executable(tridec tridec.c)
target_link_libraries(tridec utils)

install(TARGETS tridec DESTINATION bin)
