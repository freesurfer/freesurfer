project(rbftest)

include_directories(${FS_INCLUDE_DIRS})

add_executable(rbftest rbftest.c)
target_link_libraries(rbftest utils)

install(TARGETS rbftest DESTINATION bin)
