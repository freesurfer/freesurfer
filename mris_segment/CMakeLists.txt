project(mris_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_segment mris_segment.cpp)
target_link_libraries(mris_segment utils)

install(TARGETS mris_segment DESTINATION bin)
