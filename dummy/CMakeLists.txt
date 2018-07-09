project(dummy)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(dummy dummy.c)
target_link_libraries(dummy utils)

install(TARGETS dummy DESTINATION bin)
