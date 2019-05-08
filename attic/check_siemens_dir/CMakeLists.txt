project(check_siemens_dir)

include_directories(${FS_INCLUDE_DIRS})

add_executable(check_siemens_dir check_siemens_dir.cpp)
target_link_libraries(check_siemens_dir utils)

install(TARGETS check_siemens_dir DESTINATION bin)
