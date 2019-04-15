project(mrisp_paint)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mrisp_paint mrisp_paint.cpp)
add_help(mrisp_paint mrisp_paint.help.xml)
target_link_libraries(mrisp_paint utils)

install(TARGETS mrisp_paint DESTINATION bin)
