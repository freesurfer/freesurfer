project(mrisp_paint)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mrisp_paint mrisp_paint.c)
add_help(mrisp_paint mrisp_paint.help.xml)
target_link_libraries(mrisp_paint utils)

install(TARGETS mrisp_paint DESTINATION bin)
