project(tkregister2)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(tkregister2 tkregister2.c)
target_link_libraries(tkregister2 utils)

install(TARGETS tkregister2 DESTINATION bin)
