project(mris_copy_header)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_copy_header mris_copy_header.c)
target_link_libraries(mris_copy_header utils)

install(TARGETS mris_copy_header DESTINATION bin)
