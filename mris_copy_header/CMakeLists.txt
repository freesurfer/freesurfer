project(mris_copy_header)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_copy_header mris_copy_header.cpp)
target_link_libraries(mris_copy_header utils)

install(TARGETS mris_copy_header DESTINATION bin)
