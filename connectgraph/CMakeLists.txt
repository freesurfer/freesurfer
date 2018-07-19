project(connectgraph)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(connectgraph connectgraph.cpp MyCmdLineParser.cpp)
target_link_libraries(connectgraph utils)

install(TARGETS connectgraph DESTINATION bin)
