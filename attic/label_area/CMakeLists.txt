project(label_area)

include_directories(${FS_INCLUDE_DIRS})

add_executable(label_area label_area.cpp)
target_link_libraries(label_area utils)

install(TARGETS label_area DESTINATION bin)
