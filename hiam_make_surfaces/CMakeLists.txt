project(hiam_make_surfaces)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(hiam_make_surfaces hiam_make_surfaces.c)
target_link_libraries(hiam_make_surfaces utils)

install(TARGETS hiam_make_surfaces DESTINATION bin)
