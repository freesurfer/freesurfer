project(tridec)

include_directories(${FS_INCLUDE_DIRS})

add_executable(tridec tridec.cpp)
target_link_libraries(tridec utils)

install(TARGETS tridec DESTINATION bin)
