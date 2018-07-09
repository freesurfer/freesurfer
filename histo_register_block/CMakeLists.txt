project(histo_register_block)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(histo_register_block histo_register_block.c)
target_link_libraries(histo_register_block utils)

install(TARGETS histo_register_block DESTINATION bin)
