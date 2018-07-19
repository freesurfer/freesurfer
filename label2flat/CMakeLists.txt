project(label2flat)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(label2flat label2flat.c)
target_link_libraries(label2flat utils)

install(TARGETS label2flat DESTINATION bin)
