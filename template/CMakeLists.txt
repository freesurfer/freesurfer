project(template)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(template template.c)
target_link_libraries(template utils)

install(TARGETS template DESTINATION bin)
