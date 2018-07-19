project(check_siemens_dir)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(check_siemens_dir check_siemens_dir.c)
target_link_libraries(check_siemens_dir utils)

install(TARGETS check_siemens_dir DESTINATION bin)
