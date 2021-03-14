project(io)

include_directories(${FS_INCLUDE_DIRS})

add_executable(fs_temp_file fs_temp_file.cpp)
add_help(fs_temp_file fs_temp_file.help.xml)
target_link_libraries(fs_temp_file utils)

add_executable(fs_temp_dir fs_temp_dir.cpp)
add_help(fs_temp_dir fs_temp_dir.help.xml)
target_link_libraries(fs_temp_dir utils)

install(TARGETS fs_temp_file fs_temp_dir DESTINATION bin)
