project(mrisp_write)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mrisp_write mrisp_write.c)
target_link_libraries(mrisp_write utils)
install(TARGETS mrisp_write DESTINATION bin)

