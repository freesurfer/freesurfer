project(rbftest)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(rbftest rbftest.c)
target_link_libraries(rbftest utils)

install(TARGETS rbftest DESTINATION bin)
