project(tetgen)

include_directories(${FS_INCLUDE_DIRS})

add_executable(tetgen tetgen.cxx)
target_link_libraries(tetgen utils)

install(TARGETS tetgen DESTINATION bin)
