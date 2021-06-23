project(pointset2label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(pointset2label pointset2label.cpp)
target_link_libraries(pointset2label utils)

install(TARGETS pointset2label DESTINATION bin)
