project(histo_segment)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(histo_segment histo_segment.c)
target_link_libraries(histo_segment utils)

install(TARGETS histo_segment DESTINATION bin)
