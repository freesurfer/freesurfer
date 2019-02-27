project(test_makevol)

include_directories(${FS_INCLUDE_DIRS})

add_executable(makevol main.cpp)
target_link_libraries(makevol utils)

install(TARGETS makevol DESTINATION bin)