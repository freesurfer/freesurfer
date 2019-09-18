project(histo_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(histo_segment histo_segment.cpp)
target_link_libraries(histo_segment utils)

install(TARGETS histo_segment DESTINATION bin)
