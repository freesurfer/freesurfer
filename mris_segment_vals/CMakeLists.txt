project(mris_segment_vals)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_segment_vals mris_segment_vals.c)
target_link_libraries(mris_segment_vals utils)

install(TARGETS mris_segment_vals DESTINATION bin)
