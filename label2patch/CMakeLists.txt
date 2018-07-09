project(label2patch)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(label2patch label2patch.c)
target_link_libraries(label2patch utils)

install(TARGETS label2patch DESTINATION bin)
