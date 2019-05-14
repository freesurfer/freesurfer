project(rbftest)

include_directories(${FS_INCLUDE_DIRS})

add_executable(rbftest rbftest.cpp)
target_link_libraries(rbftest utils)

install(TARGETS rbftest DESTINATION bin)
