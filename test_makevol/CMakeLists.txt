project(test_makevol)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(makevol main.c)
target_link_libraries(makevol utils)

install(TARGETS makevol DESTINATION bin)