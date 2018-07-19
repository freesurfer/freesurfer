project(tetgen)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(tetgen tetgen.cxx)
target_link_libraries(tetgen utils)

install(TARGETS tetgen DESTINATION bin)
