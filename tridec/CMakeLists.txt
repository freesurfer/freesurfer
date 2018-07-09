project(tridec)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(tridec tridec.c)
target_link_libraries(tridec utils)

install(TARGETS tridec DESTINATION bin)
