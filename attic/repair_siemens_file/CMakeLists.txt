project(repair_siemens_file)

include_directories(${FS_INCLUDE_DIRS})

add_executable(repair_siemens_file repair_siemens_file.cpp)
target_link_libraries(repair_siemens_file utils)

install(TARGETS repair_siemens_file DESTINATION bin)
