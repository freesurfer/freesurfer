project(hiam_make_surfaces)

include_directories(${FS_INCLUDE_DIRS})

add_executable(hiam_make_surfaces hiam_make_surfaces.cpp)
target_link_libraries(hiam_make_surfaces utils)

install(TARGETS hiam_make_surfaces DESTINATION bin)
