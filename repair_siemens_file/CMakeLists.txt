project(repair_siemens_file)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(repair_siemens_file repair_siemens_file.c)
target_link_libraries(repair_siemens_file utils)

install(TARGETS repair_siemens_file DESTINATION bin)
