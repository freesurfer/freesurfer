project(unittests)

include_directories(${FS_INCLUDE_DIRS})

add_executable(testchklc testchklc.cpp)
target_link_libraries(testchklc utils)
