project(mris_segment)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_segment mris_segment.c)
target_link_libraries(mris_segment utils)

install(TARGETS mris_segment DESTINATION bin)
