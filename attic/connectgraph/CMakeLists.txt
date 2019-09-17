project(connectgraph)

include_directories(${FS_INCLUDE_DIRS})

add_executable(connectgraph connectgraph.cpp MyCmdLineParser.cpp)
target_link_libraries(connectgraph utils)

install(TARGETS connectgraph DESTINATION bin)
