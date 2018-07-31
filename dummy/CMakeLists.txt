project(dummy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(dummy dummy.c)
target_link_libraries(dummy utils)

install(TARGETS dummy DESTINATION bin)
