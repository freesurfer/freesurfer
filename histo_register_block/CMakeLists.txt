project(histo_register_block)

include_directories(${FS_INCLUDE_DIRS})

add_executable(histo_register_block histo_register_block.cpp)
target_link_libraries(histo_register_block utils)

install(TARGETS histo_register_block DESTINATION bin)
