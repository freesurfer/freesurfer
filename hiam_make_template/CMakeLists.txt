project(hiam_make_template)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(hiam_make_template hiam_make_template.c)
target_link_libraries(hiam_make_template utils)

install(TARGETS hiam_make_template DESTINATION bin)
