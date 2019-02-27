project(mrisp_write)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mrisp_write mrisp_write.cpp)
add_help(mrisp_write mrisp_write.help.xml)
target_link_libraries(mrisp_write utils)
install(TARGETS mrisp_write DESTINATION bin)

