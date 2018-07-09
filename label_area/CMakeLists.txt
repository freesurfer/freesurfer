project(label_area)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(label_area label_area.c)
target_link_libraries(label_area utils)

install(TARGETS label_area DESTINATION bin)
