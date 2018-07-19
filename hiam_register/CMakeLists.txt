project(hiam_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(hiam_register hiam_register.c)
target_link_libraries(hiam_register utils)

install(TARGETS hiam_register DESTINATION bin)
