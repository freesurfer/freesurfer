project(hiam_make_template)

include_directories(${FS_INCLUDE_DIRS})

add_executable(hiam_make_template hiam_make_template.cpp)
target_link_libraries(hiam_make_template utils)

install(TARGETS hiam_make_template DESTINATION bin)
