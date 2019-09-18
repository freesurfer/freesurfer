project(label2patch)

include_directories(${FS_INCLUDE_DIRS})

add_executable(label2patch label2patch.cpp)
target_link_libraries(label2patch utils)

install(TARGETS label2patch DESTINATION bin)
