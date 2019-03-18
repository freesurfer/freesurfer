project(mris_segment_vals)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_segment_vals mris_segment_vals.cpp)
target_link_libraries(mris_segment_vals utils)

install(TARGETS mris_segment_vals DESTINATION bin)
